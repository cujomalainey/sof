// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2019-2021 Intel Corporation. All rights reserved.
//
// Author: Tomasz Lauda <tomasz.lauda@linux.intel.com>

#include <sof/drivers/timer.h>
#include <sof/lib/alloc.h>
#include <sof/lib/cpu.h>
#include <sof/lib/memory.h>
#include <sof/math/numbers.h>
#include <sof/platform.h>
#include <sof/schedule/ll_schedule.h>
#include <sof/schedule/ll_schedule_domain.h>
#include <sof/schedule/schedule.h>
#include <sof/schedule/task.h>
#include <ipc/topology.h>
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <kernel.h>
#include <sys_clock.h>

/*
 * Currently the Zephyr clock rate is part it's Kconfig known at build time.
 * SOF on Intel CAVS platforms currently only aligns with Zephyr when both
 * use the CAVS 19.2 MHz SSP clock. TODO - needs runtime alignment.
 */
#if CONFIG_XTENSA && CONFIG_CAVS && !CONFIG_CAVS_TIMER
#error "Zephyr uses 19.2MHz clock derived from SSP which must be enabled."
#endif

#define ZEPHYR_LL_STACK_SIZE	8192

#define LL_TIMER_PERIOD_US 1000 /* period in microseconds */

K_KERNEL_STACK_ARRAY_DEFINE(ll_sched_stack, CONFIG_CORE_COUNT, ZEPHYR_LL_STACK_SIZE);

struct zephyr_domain_thread {
	struct k_thread ll_thread;
	void (*handler)(void *arg);
	void *arg;
};

struct zephyr_domain {
	struct k_timer timer;
	struct k_sem sem;
	struct timer *ll_timer;
	struct zephyr_domain_thread domain_thread[CONFIG_CORE_COUNT];
	struct ll_schedule_domain *ll_domain;
};

static void zephyr_domain_thread_fn(void *p1, void *p2, void *p3)
{
	struct zephyr_domain *zephyr_domain = p1;
	int core = cpu_get_id();
	struct zephyr_domain_thread *dt = zephyr_domain->domain_thread + core;

	for (;;) {
		/* immediately go to sleep, waiting to be woken up by the timer */
		k_sem_take(&zephyr_domain->sem, K_FOREVER);

		dt->handler(dt->arg);
	}
}

/* Timer callback: runs in timer IRQ context */
static void zephyr_domain_timer_fn(struct k_timer *timer)
{
	struct zephyr_domain *zephyr_domain = timer->user_data;
	int core;

	if (!zephyr_domain)
		return;

	for (core = 0; core < CONFIG_CORE_COUNT; core++)
		if (zephyr_domain->domain_thread[core].handler)
			k_sem_give(&zephyr_domain->sem);
}

static int zephyr_domain_register(struct ll_schedule_domain *domain,
				  uint64_t period, struct task *task,
				  void (*handler)(void *arg), void *arg)
{
	struct zephyr_domain *zephyr_domain = ll_sch_domain_get_pdata(domain);
	int core = cpu_get_id();
	struct zephyr_domain_thread *dt = zephyr_domain->domain_thread + core;
	char thread_name[] = "ll_thread0";
	k_tid_t thread;

	tr_dbg(&ll_tr, "zephyr_domain_register()");

	/* domain work only needs registered once on each core */
	if (dt->handler)
		return 0;

	dt->handler = handler;
	dt->arg = arg;

	thread_name[sizeof(thread_name) - 2] = '0' + core;

	thread = k_thread_create(&dt->ll_thread,
				 ll_sched_stack[core],
				 ZEPHYR_LL_STACK_SIZE,
				 zephyr_domain_thread_fn, zephyr_domain, NULL, NULL,
				 -CONFIG_NUM_COOP_PRIORITIES, 0, K_FOREVER);

	k_thread_cpu_mask_clear(thread);
	k_thread_cpu_mask_enable(thread, core);
	k_thread_name_set(thread, thread_name);

	k_thread_start(thread);

	if (!zephyr_domain->timer.user_data) {
		k_timeout_t start = {0};

		k_timer_init(&zephyr_domain->timer, zephyr_domain_timer_fn, NULL);
		zephyr_domain->timer.user_data = zephyr_domain;

		k_timer_start(&zephyr_domain->timer, start, K_USEC(LL_TIMER_PERIOD_US));
	}

	tr_info(&ll_tr, "zephyr_domain_register domain->type %d domain->clk %d domain->ticks_per_ms %d period %d",
		domain->type, domain->clk, domain->ticks_per_ms, (uint32_t)LL_TIMER_PERIOD_US);

	return 0;
}

static int zephyr_domain_unregister(struct ll_schedule_domain *domain,
				   struct task *task, uint32_t num_tasks)
{
	struct zephyr_domain *zephyr_domain = ll_sch_domain_get_pdata(domain);
	int core = cpu_get_id();

	tr_dbg(&ll_tr, "zephyr_domain_unregister()");

	/* tasks still registered on this core */
	if (num_tasks)
		return 0;

	if (!atomic_read(&domain->total_num_tasks))
		k_timer_stop(&zephyr_domain->timer);

	k_thread_abort(&zephyr_domain->domain_thread[core].ll_thread);
	zephyr_domain->domain_thread[core].handler = NULL;

	tr_info(&ll_tr, "zephyr_domain_unregister domain->type %d domain->clk %d",
		domain->type, domain->clk);

	return 0;
}

static bool zephyr_domain_is_pending(struct ll_schedule_domain *domain,
				     struct task *task, struct comp_dev **comp)
{
	struct zephyr_domain *zephyr_domain = ll_sch_domain_get_pdata(domain);

	return task->start <= platform_timer_get_atomic(zephyr_domain->ll_timer);
}

static const struct ll_schedule_domain_ops zephyr_domain_ops = {
	.domain_register	= zephyr_domain_register,
	.domain_unregister	= zephyr_domain_unregister,
	.domain_is_pending	= zephyr_domain_is_pending
};

struct ll_schedule_domain *timer_domain_init(struct timer *timer, int clk)
{
	struct ll_schedule_domain *domain;
	struct zephyr_domain *zephyr_domain;

	domain = domain_init(SOF_SCHEDULE_LL_TIMER, clk, false,
			     &zephyr_domain_ops);

	zephyr_domain = rzalloc(SOF_MEM_ZONE_SYS_SHARED, 0, SOF_MEM_CAPS_RAM,
				sizeof(*zephyr_domain));

	zephyr_domain->ll_timer = timer;
	zephyr_domain->ll_domain = domain;
	/* 10 is rather random, we better not accumulate 10 missed timer interrupts */
	k_sem_init(&zephyr_domain->sem, 0, 10);

	ll_sch_domain_set_pdata(domain, zephyr_domain);

	return domain;
}
