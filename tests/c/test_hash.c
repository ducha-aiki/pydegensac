#include <assert.h>
#include "hash.h"

int main(void)
{
    HashTable table, other;
    unsigned i, round;
    const unsigned count = 3 * HT_BLOCK_CAPACITY + 7;

    htInit(&table);
    htInit(&other);
    htInsert(&other, 5 * HT_FIELDS, 4, 999);
    for (round = 0; round < 3; ++round) {
        /* Collisions and enough entries to cross multiple allocation blocks. */
        for (i = 0; i < count; ++i)
            htInsert(&table, (i + 1) * HT_FIELDS, i % 17, i);
        for (i = 0; i < count; ++i) {
            assert(htContains(&table, (i + 1) * HT_FIELDS, i % 17, i) == (int)i);
            assert(htContains(&table, (i + 1) * HT_FIELDS, i % 17, -2) == (int)i);
        }
        htInsert(&table, 5 * HT_FIELDS, 4, 900);
        htInsert(&table, 5 * HT_FIELDS, 4, 901);
        /* Prefer the requested iteration, otherwise the most recent match. */
        assert(htContains(&table, 5 * HT_FIELDS, 4, 4) == 4);
        assert(htContains(&table, 5 * HT_FIELDS, 4, 900) == 900);
        assert(htContains(&table, 5 * HT_FIELDS, 4, -2) == 901);
        assert(htContains(&table, 5 * HT_FIELDS, 999, 901) == -1);
        htClear(&table);
        htClear(&table);
        for (i = 0; i < count; ++i)
            assert(htContains(&table, (i + 1) * HT_FIELDS, i % 17, i) == -1);
        assert(htContains(&other, 5 * HT_FIELDS, 4, 999) == 999);
    }
    htClear(&other);
    return 0;
}
