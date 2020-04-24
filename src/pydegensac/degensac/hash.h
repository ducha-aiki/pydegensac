#include <stdint.h>  /* Replace with <stdint.h if appropriate */
#include <stdlib.h>

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

#define HT_FIELDS 64

uint32_t SuperFastHash (const char * data, int len);

typedef struct _HashField
{
    uint32_t hash;
    int length;
    int iterID;
    double thr;
    struct _HashField * next;
} HashField;

typedef struct _HashTable
{
    HashField * fields[HT_FIELDS];
} HashTable;
HashTable HASH_TABLE;

void htInit(HashTable * ht);

void htClear(HashTable * ht);

void htInsert(HashTable * ht, uint32_t hash, int length, int iterID);

int htContains(HashTable * ht, uint32_t hash, int length, int iterID);
