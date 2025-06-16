/* SPDX-License-Identifier: MIT */
#ifndef LIBURING_COMPAT_H
#define LIBURING_COMPAT_H

#include <stdint.h>

struct __kernel_timespec {
	int64_t		tv_sec;
	long long	tv_nsec;
};

/* <linux/time_types.h> is not available, so it can't be included */
#define UAPI_LINUX_IO_URING_H_SKIP_LINUX_TIME_TYPES_H 1

#include <inttypes.h>

struct open_how {
	uint64_t	flags;
	uint64_t	mode;
	uint64_t	resolve;
};

#endif
