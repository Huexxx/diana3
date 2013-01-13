/*
 * mm/readahead.c - address_space-level file readahead.
 *
 * Copyright (C) 2002, Linus Torvalds
 *
 * 09Apr2002	Andrew Morton
 *		Initial version.
 */

/*
 * Notes on readahead size.
 *
 * The default max readahead size is VM_MAX_READAHEAD=512k,
 * which can be changed by user with boot time parameter "readahead="
 * or runtime interface "/sys/devices/virtual/bdi/default/read_ahead_kb".
 * The latter normally only takes effect in future for hot added devices.
 *
 * The effective max readahead size for each block device can be accessed with
 * 1) the `blockdev` command
 * 2) /sys/block/sda/queue/read_ahead_kb
 * 3) /sys/devices/virtual/bdi/$(env stat -c '%t:%T' /dev/sda)/read_ahead_kb
 *
 * They are typically initialized with the global default size, however may be
 * auto scaled down for small devices in add_disk(). NFS, software RAID, btrfs
 * etc. have special rules to setup their default readahead size.
 *
 * The mmap read-around size typically equals with readahead size, with an
 * extra limit proportional to system memory size.  For example, a 64MB box
 * will have a 64KB read-around size limit, 128MB mem => 128KB limit, etc.
 */

#include <linux/kernel.h>
#include <linux/fs.h>
#include <linux/memcontrol.h>
#include <linux/gfp.h>
#include <linux/mm.h>
#include <linux/mm_inline.h>
#include <linux/module.h>
#include <linux/blkdev.h>
#include <linux/backing-dev.h>
#include <linux/task_io_accounting_ops.h>
#include <linux/pagevec.h>
#include <linux/pagemap.h>

#define CREATE_TRACE_POINTS
#include <trace/events/readahead.h>

/*
 * Set async size to 1/# of the thrashing threshold.
 */
#define READAHEAD_ASYNC_RATIO	8

#define MIN_READAHEAD_PAGES DIV_ROUND_UP(VM_MIN_READAHEAD*1024, PAGE_CACHE_SIZE)

static int __init config_readahead_size(char *str)
{
	unsigned long bytes;

	if (!str)
		return -EINVAL;
	bytes = memparse(str, &str);
	if (*str != '\0')
		return -EINVAL;

	if (bytes) {
		if (bytes < PAGE_CACHE_SIZE)	/* missed 'k'/'m' suffixes? */
			return -EINVAL;
		if (bytes > 128 << 20)		/* limit to 128MB */
			bytes = 128 << 20;
	}

	default_backing_dev_info.ra_pages = bytes / PAGE_CACHE_SIZE;
	return 0;
}
early_param("readahead", config_readahead_size);

/*
 * Initialise a struct file's readahead state.  Assumes that the caller has
 * memset *ra to zero.
 */
void
file_ra_state_init(struct file_ra_state *ra, struct address_space *mapping)
{
	ra->ra_pages = mapping->backing_dev_info->ra_pages;
	ra->prev_pos = -1;
}
EXPORT_SYMBOL_GPL(file_ra_state_init);

#define list_to_page(head) (list_entry((head)->prev, struct page, lru))

#ifdef CONFIG_READAHEAD_STATS
#include <linux/seq_file.h>
#include <linux/debugfs.h>
enum ra_account {
	/* number of readaheads */
	RA_ACCOUNT_COUNT,	/* readahead request */
	RA_ACCOUNT_EOF,		/* readahead request contains/beyond EOF page */
	RA_ACCOUNT_CHIT,	/* readahead request covers some cached pages */
	RA_ACCOUNT_IOCOUNT,	/* readahead IO */
	RA_ACCOUNT_SYNC,	/* readahead IO that is synchronous */
	RA_ACCOUNT_MMAP,	/* readahead IO by mmap accesses */
	/* number of readahead pages */
	RA_ACCOUNT_SIZE,	/* readahead size */
	RA_ACCOUNT_ASIZE,	/* readahead async size */
	RA_ACCOUNT_ACTUAL,	/* readahead actual IO size */
	/* end mark */
	RA_ACCOUNT_MAX,
};

static unsigned long ra_stats[RA_PATTERN_MAX][RA_ACCOUNT_MAX];

static void readahead_stats(struct address_space *mapping,
			    pgoff_t offset,
			    unsigned long req_size,
			    unsigned int ra_flags,
			    pgoff_t start,
			    unsigned int size,
			    unsigned int async_size,
			    int actual)
{
	unsigned int pattern = ra_pattern(ra_flags);

	ra_stats[pattern][RA_ACCOUNT_COUNT]++;
	ra_stats[pattern][RA_ACCOUNT_SIZE] += size;
	ra_stats[pattern][RA_ACCOUNT_ASIZE] += async_size;
	ra_stats[pattern][RA_ACCOUNT_ACTUAL] += actual;

	if (actual < size) {
		if (start + size >
		    (i_size_read(mapping->host) - 1) >> PAGE_CACHE_SHIFT)
			ra_stats[pattern][RA_ACCOUNT_EOF]++;
		else
			ra_stats[pattern][RA_ACCOUNT_CHIT]++;
	}

	if (!actual)
		return;

	ra_stats[pattern][RA_ACCOUNT_IOCOUNT]++;

	if (start <= offset && start + size > offset)
		ra_stats[pattern][RA_ACCOUNT_SYNC]++;

	if (ra_flags & READAHEAD_MMAP)
		ra_stats[pattern][RA_ACCOUNT_MMAP]++;
}

static int readahead_stats_show(struct seq_file *s, void *_)
{
	static const char * const ra_pattern_names[] = {
		[RA_PATTERN_INITIAL]		= "initial",
		[RA_PATTERN_SUBSEQUENT]		= "subsequent",
		[RA_PATTERN_CONTEXT]		= "context",
		[RA_PATTERN_THRASH]		= "thrash",
		[RA_PATTERN_MMAP_AROUND]	= "around",
		[RA_PATTERN_FADVISE]		= "fadvise",
		[RA_PATTERN_RANDOM]		= "random",
		[RA_PATTERN_ALL]		= "all",
	};
	unsigned long count, iocount;
	unsigned long i;

	seq_printf(s, "%-10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
			"pattern",
			"readahead", "eof_hit", "cache_hit",
			"io", "sync_io", "mmap_io",
			"size", "async_size", "io_size");

	for (i = 0; i < RA_PATTERN_MAX; i++) {
		count = ra_stats[i][RA_ACCOUNT_COUNT];
		iocount = ra_stats[i][RA_ACCOUNT_IOCOUNT];
		/*
		 * avoid division-by-zero
		 */
		if (count == 0)
			count = 1;
		if (iocount == 0)
			iocount = 1;

		seq_printf(s, "%-10s %10lu %10lu %10lu %10lu %10lu %10lu "
			   "%10lu %10lu %10lu\n",
				ra_pattern_names[i],
				ra_stats[i][RA_ACCOUNT_COUNT],
				ra_stats[i][RA_ACCOUNT_EOF],
				ra_stats[i][RA_ACCOUNT_CHIT],
				ra_stats[i][RA_ACCOUNT_IOCOUNT],
				ra_stats[i][RA_ACCOUNT_SYNC],
				ra_stats[i][RA_ACCOUNT_MMAP],
				ra_stats[i][RA_ACCOUNT_SIZE]   / count,
				ra_stats[i][RA_ACCOUNT_ASIZE]  / count,
				ra_stats[i][RA_ACCOUNT_ACTUAL] / iocount);
	}

	return 0;
}

static int readahead_stats_open(struct inode *inode, struct file *file)
{
	return single_open(file, readahead_stats_show, NULL);
}

static ssize_t readahead_stats_write(struct file *file, const char __user *buf,
				     size_t size, loff_t *offset)
{
	memset(ra_stats, 0, sizeof(ra_stats));
	return size;
}

static struct file_operations readahead_stats_fops = {
	.owner		= THIS_MODULE,
	.open		= readahead_stats_open,
	.write		= readahead_stats_write,
	.read		= seq_read,
	.llseek		= seq_lseek,
	.release	= single_release,
};

static struct dentry *ra_debug_root;

static int debugfs_create_readahead(void)
{
	struct dentry *debugfs_stats;

	ra_debug_root = debugfs_create_dir("readahead", NULL);
	if (!ra_debug_root)
		goto out;

	debugfs_stats = debugfs_create_file("stats", 0644, ra_debug_root,
					    NULL, &readahead_stats_fops);
	if (!debugfs_stats)
		goto out;

	return 0;
out:
	printk(KERN_ERR "readahead: failed to create debugfs entries\n");
	return -ENOMEM;
}

static int __init readahead_init(void)
{
	debugfs_create_readahead();
	return 0;
}

static void __exit readahead_exit(void)
{
	debugfs_remove_recursive(ra_debug_root);
}

module_init(readahead_init);
module_exit(readahead_exit);
#endif

static void readahead_event(struct address_space *mapping,
			    pgoff_t offset,
			    unsigned long req_size,
			    unsigned int ra_flags,
			    pgoff_t start,
			    unsigned int size,
			    unsigned int async_size,
			    unsigned int actual)
{
#ifdef CONFIG_READAHEAD_STATS
	readahead_stats(mapping, offset, req_size, ra_flags,
			start, size, async_size, actual);
	readahead_stats(mapping, offset, req_size,
			RA_PATTERN_ALL << READAHEAD_PATTERN_SHIFT,
			start, size, async_size, actual);
#endif
	trace_readahead(mapping, offset, req_size, ra_flags,
			start, size, async_size, actual);
}

/*
 * see if a page needs releasing upon read_cache_pages() failure
 * - the caller of read_cache_pages() may have set PG_private or PG_fscache
 *   before calling, such as the NFS fs marking pages that are cached locally
 *   on disk, thus we need to give the fs a chance to clean up in the event of
 *   an error
 */
static void read_cache_pages_invalidate_page(struct address_space *mapping,
					     struct page *page)
{
	if (page_has_private(page)) {
		if (!trylock_page(page))
			BUG();
		page->mapping = mapping;
		do_invalidatepage(page, 0);
		page->mapping = NULL;
		unlock_page(page);
	}
	page_cache_release(page);
}

/*
 * release a list of pages, invalidating them first if need be
 */
static void read_cache_pages_invalidate_pages(struct address_space *mapping,
					      struct list_head *pages)
{
	struct page *victim;

	while (!list_empty(pages)) {
		victim = list_to_page(pages);
		list_del(&victim->lru);
		read_cache_pages_invalidate_page(mapping, victim);
	}
}

/**
 * read_cache_pages - populate an address space with some pages & start reads against them
 * @mapping: the address_space
 * @pages: The address of a list_head which contains the target pages.  These
 *   pages have their ->index populated and are otherwise uninitialised.
 * @filler: callback routine for filling a single page.
 * @data: private data for the callback routine.
 *
 * Hides the details of the LRU cache etc from the filesystems.
 */
int read_cache_pages(struct address_space *mapping, struct list_head *pages,
			int (*filler)(void *, struct page *), void *data)
{
	struct page *page;
	int ret = 0;

	while (!list_empty(pages)) {
		page = list_to_page(pages);
		list_del(&page->lru);
		if (add_to_page_cache_lru(page, mapping,
					page->index, GFP_KERNEL)) {
			read_cache_pages_invalidate_page(mapping, page);
			continue;
		}
		page_cache_release(page);

		ret = filler(data, page);
		if (unlikely(ret)) {
			read_cache_pages_invalidate_pages(mapping, pages);
			break;
		}
		task_io_account_read(PAGE_CACHE_SIZE);
	}
	return ret;
}

EXPORT_SYMBOL(read_cache_pages);

static int read_pages(struct address_space *mapping, struct file *filp,
		struct list_head *pages, unsigned nr_pages)
{
	struct blk_plug plug;
	unsigned page_idx;
	int ret;

	blk_start_plug(&plug);

	if (mapping->a_ops->readpages) {
		ret = mapping->a_ops->readpages(filp, mapping, pages, nr_pages);
		/* Clean up the remaining pages */
		put_pages_list(pages);
		goto out;
	}

	for (page_idx = 0; page_idx < nr_pages; page_idx++) {
		struct page *page = list_to_page(pages);
		list_del(&page->lru);
		if (!add_to_page_cache_lru(page, mapping,
					page->index, GFP_KERNEL)) {
			mapping->a_ops->readpage(filp, page);
		}
		page_cache_release(page);
	}
	ret = 0;

out:
	blk_finish_plug(&plug);

	return ret;
}

/*
 * The file range is expected to be accessed in near future.  Move pages
 * (possibly in inactive lru tail) to lru head, so that they are retained
 * in memory for some reasonable time.
 */
static void retain_inactive_pages(struct address_space *mapping,
				  pgoff_t index, int len)
{
	int i;
	struct page *page;
	struct zone *zone;

	for (i = 0; i < len; i++) {
		page = find_get_page(mapping, index + i);
		if (!page)
			continue;

		zone = page_zone(page);
		spin_lock_irq(&zone->lru_lock);

		if (PageLRU(page) &&
		    !PageActive(page) &&
		    !PageUnevictable(page)) {
			int lru = page_lru_base_type(page);

			del_page_from_lru_list(zone, page, lru);
			add_page_to_lru_list(zone, page, lru);
		}

		spin_unlock_irq(&zone->lru_lock);
		put_page(page);
	}
}

/*
 * __do_page_cache_readahead() actually reads a chunk of disk.  It allocates all
 * the pages first, then submits them all for I/O. This avoids the very bad
 * behaviour which would occur if page allocations are causing VM writeback.
 * We really don't want to intermingle reads and writes like that.
 *
 * Returns the number of pages requested, or the maximum amount of I/O allowed.
 */
static int
__do_page_cache_readahead(struct address_space *mapping, struct file *filp,
			pgoff_t offset, unsigned long nr_to_read,
			unsigned long lookahead_size)
{
	struct inode *inode = mapping->host;
	struct page *page;
	unsigned long end_index;	/* The last page we want to read */
	LIST_HEAD(page_pool);
	int page_idx;
	int ret = 0;
	loff_t isize = i_size_read(inode);

	if (isize == 0)
		goto out;

	end_index = ((isize - 1) >> PAGE_CACHE_SHIFT);

	/*
	 * Preallocate as many pages as we will need.
	 */
	for (page_idx = 0; page_idx < nr_to_read; page_idx++) {
		pgoff_t page_offset = offset + page_idx;

		if (page_offset > end_index)
			break;

		rcu_read_lock();
		page = radix_tree_lookup(&mapping->page_tree, page_offset);
		rcu_read_unlock();
		if (page)
			continue;

		page = page_cache_alloc_readahead(mapping);
		if (!page)
			break;
		page->index = page_offset;
		list_add(&page->lru, &page_pool);
		if (page_idx == nr_to_read - lookahead_size)
			SetPageReadahead(page);
		ret++;
	}

	/*
	 * Normally readahead will auto stop on cached segments, so we won't
	 * hit many cached pages. If it does happen, bring the inactive pages
	 * adjecent to the newly prefetched ones(if any).
	 */
	if (ret < nr_to_read)
		retain_inactive_pages(mapping, offset, page_idx);

	/*
	 * Now start the IO.  We ignore I/O errors - if the page is not
	 * uptodate then the caller will launch readpage again, and
	 * will then handle the error.
	 */
	if (ret)
		read_pages(mapping, filp, &page_pool, ret);
	BUG_ON(!list_empty(&page_pool));
out:
	return ret;
}

/*
 * Chunk the readahead into 2 megabyte units, so that we don't pin too much
 * memory at once.
 */
int force_page_cache_readahead(struct address_space *mapping, struct file *filp,
		pgoff_t offset, unsigned long nr_to_read)
{
	int ret = 0;

	if (unlikely(!mapping->a_ops->readpage && !mapping->a_ops->readpages))
		return -EINVAL;

	nr_to_read = max_sane_readahead(nr_to_read);
	while (nr_to_read) {
		int err;

		unsigned long this_chunk = (2 * 1024 * 1024) / PAGE_CACHE_SIZE;

		if (this_chunk > nr_to_read)
			this_chunk = nr_to_read;
		err = __do_page_cache_readahead(mapping, filp,
						offset, this_chunk, 0);
		if (err < 0) {
			ret = err;
			break;
		}
		ret += err;
		offset += this_chunk;
		nr_to_read -= this_chunk;
	}

	readahead_event(mapping, offset, nr_to_read,
			RA_PATTERN_FADVISE << READAHEAD_PATTERN_SHIFT,
			offset, nr_to_read, 0, ret);

	return ret;
}

/*
 * Given a desired number of PAGE_CACHE_SIZE readahead pages, return a
 * sensible upper limit.
 */
unsigned long max_sane_readahead(unsigned long nr)
{
	return min(nr, (node_page_state(numa_node_id(), NR_INACTIVE_FILE)
		+ node_page_state(numa_node_id(), NR_FREE_PAGES)) / 2);
}

/*
 * Submit IO for the read-ahead request in file_ra_state.
 */
unsigned long ra_submit(struct file_ra_state *ra,
			struct address_space *mapping,
			struct file *filp,
			pgoff_t offset,
			unsigned long req_size)
{
	int actual;

	actual = __do_page_cache_readahead(mapping, filp,
					ra->start, ra->size, ra->async_size);

	readahead_event(mapping, offset, req_size, ra->ra_flags,
			ra->start, ra->size, ra->async_size, actual);

	return actual;
}

/*
 * Set the initial window size, round to next power of 2 and square
 * for small size, x 4 for medium, and x 2 for large
 * for 128k (32 page) max ra
 * 1-8 page = 32k initial, > 8 page = 128k initial
 */
static unsigned long get_init_ra_size(unsigned long size, unsigned long max)
{
	unsigned long newsize = roundup_pow_of_two(size);

	if (newsize <= max / 32)
		newsize = newsize * 4;
	else if (newsize <= max / 4)
		newsize = newsize * 2;
	else
		newsize = max;

	return newsize;
}

/*
 *  Get the previous window size, ramp it up, and
 *  return it as the new window size.
 */
static unsigned long get_next_ra_size(struct file_ra_state *ra,
						unsigned long max)
{
	unsigned long cur = ra->size;
	unsigned long newsize;

	if (cur < max / 16)
		newsize = 4 * cur;
	else
		newsize = 2 * cur;

	return min(newsize, max);
}

/*
 * On-demand readahead design.
 *
 * The fields in struct file_ra_state represent the most-recently-executed
 * readahead attempt:
 *
 *                        |<----- async_size ---------|
 *     |------------------- size -------------------->|
 *     |==================#===========================|
 *     ^start             ^page marked with PG_readahead
 *
 * To overlap application thinking time and disk I/O time, we do
 * `readahead pipelining': Do not wait until the application consumed all
 * readahead pages and stalled on the missing page at readahead_index;
 * Instead, submit an asynchronous readahead I/O as soon as there are
 * only async_size pages left in the readahead window. Normally async_size
 * will be equal to size, for maximum pipelining.
 *
 * In interleaved sequential reads, concurrent streams on the same fd can
 * be invalidating each other's readahead state. So we flag the new readahead
 * page at (start+size-async_size) with PG_readahead, and use it as readahead
 * indicator. The flag won't be set on already cached pages, to avoid the
 * readahead-for-nothing fuss, saving pointless page cache lookups.
 *
 * prev_pos tracks the last visited byte in the _previous_ read request.
 * It should be maintained by the caller, and will be used for detecting
 * small random reads. Note that the readahead algorithm checks loosely
 * for sequential patterns. Hence interleaved reads might be served as
 * sequential ones.
 *
 * There is a special-case: if the first page which the application tries to
 * read happens to be the first page of the file, it is assumed that a linear
 * read is about to happen and the window is immediately set to the initial size
 * based on I/O request size and the max_readahead.
 *
 * The code ramps up the readahead size aggressively at first, but slow down as
 * it approaches max_readhead.
 */

/*
 * Count contiguously cached pages from @offset-1 to @offset-@max,
 * this count is a conservative estimation of
 * 	- length of the sequential read sequence, or
 * 	- thrashing threshold in memory tight systems
 */
static pgoff_t count_history_pages(struct address_space *mapping,
				   struct file_ra_state *ra,
				   pgoff_t offset, unsigned long max)
{
	pgoff_t head;

	rcu_read_lock();
	head = radix_tree_prev_hole(&mapping->page_tree, offset - 1, max);
	rcu_read_unlock();

	return offset - 1 - head;
}

/*
 * Is @index recently readahead but not yet read by application?
 * The low boundary is permissively estimated.
 */
static bool ra_thrashed(struct file_ra_state *ra, pgoff_t index)
{
	return (index >= ra->start - ra->size &&
		index <  ra->start + ra->size);
}


/*
 * A minimal readahead algorithm for trivial sequential/random reads.
 */
static unsigned long
ondemand_readahead(struct address_space *mapping,
		   struct file_ra_state *ra, struct file *filp,
		   bool hit_readahead_marker, pgoff_t offset,
		   unsigned long req_size)
{
	unsigned long max = max_sane_readahead(ra->ra_pages);
	unsigned long tt;  /* thrashing shreshold */
	pgoff_t start;

	/*
	 * start of file
	 */
	if (!offset) {
		ra_set_pattern(ra, RA_PATTERN_INITIAL);
		ra->start = offset;
		if ((ra->ra_flags & READAHEAD_LSEEK) && req_size <= max) {
			ra->size = req_size;
			ra->async_size = 0;
			goto readit;
		}
		ra->size = get_init_ra_size(req_size, max);
		ra->async_size = ra->size > req_size ?
				 ra->size - req_size : ra->size;
		goto readit;
	}

	/*
	 * Context readahead is thrashing safe, and can adapt to near the
	 * thrashing threshold given a stable workload.
	 */
	if (ra->ra_flags & READAHEAD_THRASHED)
		goto context_readahead;

	/*
	 * It's the expected callback offset, assume sequential access.
	 * Ramp up sizes, and push forward the readahead window.
	 */
	if ((offset == (ra->start + ra->size - ra->async_size) ||
	     offset == (ra->start + ra->size))) {
		ra_set_pattern(ra, RA_PATTERN_SUBSEQUENT);
		ra->start += ra->size;
		ra->size = get_next_ra_size(ra, max);
		ra->async_size = ra->size;
		goto readit;
	}

	/*
	 * oversize read, no need to query page cache
	 */
	if (req_size > max && !hit_readahead_marker) {
		ra_set_pattern(ra, RA_PATTERN_INITIAL);
		ra->start = offset;
		ra->size = max;
		ra->async_size = max;
		goto readit;
	}

	/*
	 * page cache context based read-ahead
	 *
	 *     ==========================_____________..............
	 *                          [ current window ]
	 *                               ^offset
	 * 1)                            |---- A ---->[start
	 * 2) |<----------- H -----------|
	 * 3)                            |----------- H ----------->]end
	 *                                            [ new window ]
	 *    [=] cached,visited [_] cached,to-be-visited [.] not cached
	 *
	 * 1) A = pages ahead = previous async_size
	 * 2) H = history pages = thrashing safe size
	 * 3) H - A = new readahead size
	 */
context_readahead:
	if (hit_readahead_marker) {
		rcu_read_lock();
		start = radix_tree_next_hole(&mapping->page_tree,
					     offset + 1, max);
		rcu_read_unlock();
		/*
		 * there are enough pages ahead: no readahead
		 */
		if (!start || start - offset > max)
			return 0;
	} else
		start = offset;

	tt = count_history_pages(mapping, ra, offset,
				 READAHEAD_ASYNC_RATIO * max);
	/*
	 * no history pages cached, could be
	 * 	- a random read
	 * 	- a thrashed sequential read
	 */
	if (!tt && !hit_readahead_marker) {
		if (!ra_thrashed(ra, offset)) {
			ra_set_pattern(ra, RA_PATTERN_RANDOM);
			ra->size = min(req_size, max);
		} else {
			ra_set_pattern(ra, RA_PATTERN_THRASH);
			retain_inactive_pages(mapping, offset, min(2 * max,
						ra->start + ra->size - offset));
			ra->size = max_t(int, ra->size/2, MIN_READAHEAD_PAGES);
			ra->ra_flags |= READAHEAD_THRASHED;
		}
		ra->async_size = 0;
		ra->start = start;
		goto readit;
	}
	/*
	 * history pages start from beginning of file:
	 * it is a strong indication of long-run stream (or whole-file reads)
	 */
	if (tt >= offset)
		tt *= 2;
	/*
	 * Pages to readahead are already cached?
	 */
	if (tt <= start - offset)
		return 0;

	ra_set_pattern(ra, RA_PATTERN_CONTEXT);
	ra->start = start;
	ra->size = clamp_t(unsigned int, tt - (start - offset),
			   MIN_READAHEAD_PAGES, max);
	ra->async_size = min_t(unsigned int, ra->size,
			       1 + tt / READAHEAD_ASYNC_RATIO);

readit:
	/*
	 * Will this read hit the readahead marker made by itself?
	 * If so, trigger the readahead marker hit now, and merge
	 * the resulted next readahead window into the current one.
	 */
	if (offset == ra->start && ra->size == ra->async_size) {
		ra->async_size = get_next_ra_size(ra, max);
		ra->size += ra->async_size;
	}

	return ra_submit(ra, mapping, filp, offset, req_size);
}

/**
 * page_cache_sync_readahead - generic file readahead
 * @mapping: address_space which holds the pagecache and I/O vectors
 * @ra: file_ra_state which holds the readahead state
 * @filp: passed on to ->readpage() and ->readpages()
 * @offset: start offset into @mapping, in pagecache page-sized units
 * @req_size: hint: total size of the read which the caller is performing in
 *            pagecache pages
 *
 * page_cache_sync_readahead() should be called when a cache miss happened:
 * it will submit the read.  The readahead logic may decide to piggyback more
 * pages onto the read request if access patterns suggest it will improve
 * performance.
 */
void page_cache_sync_readahead(struct address_space *mapping,
			       struct file_ra_state *ra, struct file *filp,
			       pgoff_t offset, unsigned long req_size)
{
	/* no read-ahead */
	if (!ra->ra_pages)
		return;

	/* be dumb */
	if (filp && (filp->f_mode & FMODE_RANDOM)) {
		force_page_cache_readahead(mapping, filp, offset, req_size);
		return;
	}

	/* do read-ahead */
	ondemand_readahead(mapping, ra, filp, false, offset, req_size);
}
EXPORT_SYMBOL_GPL(page_cache_sync_readahead);

/**
 * page_cache_async_readahead - file readahead for marked pages
 * @mapping: address_space which holds the pagecache and I/O vectors
 * @ra: file_ra_state which holds the readahead state
 * @filp: passed on to ->readpage() and ->readpages()
 * @page: the page at @offset which has the PG_readahead flag set
 * @offset: start offset into @mapping, in pagecache page-sized units
 * @req_size: hint: total size of the read which the caller is performing in
 *            pagecache pages
 *
 * page_cache_async_readahead() should be called when a page is used which
 * has the PG_readahead flag; this is a marker to suggest that the application
 * has used up enough of the readahead window that we should start pulling in
 * more pages.
 */
void
page_cache_async_readahead(struct address_space *mapping,
			   struct file_ra_state *ra, struct file *filp,
			   struct page *page, pgoff_t offset,
			   unsigned long req_size)
{
	/* no read-ahead */
	if (!ra->ra_pages)
		return;

	/*
	 * Same bit is used for PG_readahead and PG_reclaim.
	 */
	if (PageWriteback(page))
		return;

	ClearPageReadahead(page);

	/*
	 * Defer asynchronous read-ahead on IO congestion.
	 */
	if (bdi_read_congested(mapping->backing_dev_info))
		return;

	/* do read-ahead */
	ondemand_readahead(mapping, ra, filp, true, offset, req_size);
}
EXPORT_SYMBOL_GPL(page_cache_async_readahead);
