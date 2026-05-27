// SPDX-License-Identifier: GPL-2.0-or-later
/**
 * @file
 * @brief Example: prints mbelib-neo version using the public API.
 */

#include <stdio.h>
#include "mbelib-neo/mbelib.h"

/**
 * @brief Program entry printing the library version to stdout.
 */
int
main(void) {
    printf("mbelib version: %s\n", mbe_versionString());
    return 0;
}
