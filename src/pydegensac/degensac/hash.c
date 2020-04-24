#include "hash.h"
//#include <mex.h>

uint32_t SuperFastHash (const char * data, int len) {
	uint32_t hash = len, tmp;
	int rem;

    if (len <= 0 || data == 0) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= ((signed char)data[sizeof (uint16_t)]) << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += (signed char)*data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

void htInit(HashTable * ht) {
	unsigned i;
	for (i = 0; i < HT_FIELDS; ++i) {
		ht->fields[i] = 0;
	}
}

void htClear(HashTable * ht) {
	unsigned i;
	HashField * hf;
	for (i = 0; i < HT_FIELDS; ++i) {
		while (ht->fields[i]) {
			hf = ht->fields[i];
			ht->fields[i] = ht->fields[i]->next;
			free(hf);
		}
	}
}

void htInsert(HashTable * ht, uint32_t hash, int length, int iterID) {
	HashField * hf = (HashField *)malloc(sizeof(HashField));
	hf->next = ht->fields[hash % HT_FIELDS];
	hf->hash = hash;
	hf->length = length;
	hf->iterID = iterID;
	ht->fields[hash % HT_FIELDS] = hf;
}

int htContains(HashTable * ht, uint32_t hash, int length, int iterID) { /* compare thrs to be sure, but slower */
	HashField * hf = ht->fields[hash % HT_FIELDS];
	while (hf) { // the same iterID
		if (hf->hash == hash && hf->length == length && hf->iterID == iterID) {
			return iterID;
		} else {
			hf = hf->next;
		}
	}
	
	hf = ht->fields[hash % HT_FIELDS];
	while (hf) { // other iterID
		if (hf->hash == hash && hf->length == length) {
			return hf->iterID;
		} else {
			hf = hf->next;
		}
	}
	return -1;
}



