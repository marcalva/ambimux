
#ifndef G_LIST_H
#define G_LIST_H

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include "kavl.h"

#define ml_declare_i(SCOPE, name, ntype, n_cmp) \
    /* struct node */ \
    typedef struct __ml_node_##name { \
        ntype data; \
        struct __ml_node_##name *next; \
    } __ml_node_##name; \
    typedef struct __ml_node_##name ml_node_##name##_t; \
    \
    /* initialize node */ \
    SCOPE __ml_node_##name *__ml_node_init_##name(void){ \
        __ml_node_##name *n = (__ml_node_##name *)calloc(1, sizeof(__ml_node_##name)); \
        if (n == NULL){ \
            fprintf(stderr, "error: __ml_node_init_%s: %s\n", #name, strerror(errno)); \
            return(NULL); \
        } \
        n->next = NULL; \
        return(n); \
    } \
    \
    /* destroy node */ \
    SCOPE __ml_node_##name *__ml_node_dstry_##name(__ml_node_##name *p){ \
        if (p == NULL) return(NULL); \
        __ml_node_##name *n = p->next; \
        free(p); \
        return n; \
    } \
    \
    /* struct list */ \
    typedef struct __ml_l_##name { \
        __ml_node_##name *head; \
        uint32_t size; \
    } ml_##name##_t; \
    \
    /* alloc list */ \
    SCOPE ml_##name##_t *ml_alloc_##name(void){ \
        ml_##name##_t *ll = (ml_##name##_t *)calloc(1, sizeof(ml_##name##_t)); \
        if (ll == NULL){ \
            fprintf(stderr, "error: ml_alloc_%s: %s\n", #name, strerror(errno)); \
            return(NULL); \
        } \
        ll->head = NULL; \
        ll->size = 0; \
        return ll; \
    } \
    \
    /* initialize list */ \
    SCOPE int ml_init_##name(ml_##name##_t *ll){ \
        if (ll == NULL){ \
            fprintf(stderr, "error: ml_init_%s: argument list is null\n", #name); \
            return(-1); \
        } \
        \
        ll->head = NULL; \
        ll->size = 0; \
        return(0); \
    }\
    \
    /* free list */ \
    SCOPE void ml_free_##name(ml_##name##_t *ll){ \
        if (ll == NULL) return; \
        __ml_node_##name *p = ll->head; \
        while (p != NULL) \
            p = __ml_node_dstry_##name(p); \
        ll->head = NULL; \
        ll->size = 0; \
    } \
    \
    /* destroy list */ \
    SCOPE void ml_dstry_##name(ml_##name##_t *ll){ \
        if (ll == NULL) return; \
        __ml_node_##name *p = ll->head; \
        while (p != NULL) \
            p = __ml_node_dstry_##name(p); \
        free(ll); \
    } \
    \
    SCOPE int ml_cpy_##name(ml_##name##_t *dest, const ml_##name##_t *src){ \
        if (dest == NULL){\
            fprintf(stderr, "error: ml_cpy_%s: dest is null\n", #name); \
            return(-1); \
        } \
        if (src == NULL){\
            fprintf(stderr, "error: ml_cpy_%s: src is null\n", #name); \
            return(-1); \
        } \
        dest->size = src->size; \
        dest->head = NULL; \
        __ml_node_##name *p_src = src->head; \
        __ml_node_##name *p_prev = NULL; \
        uint32_t n_e = 0; \
        while (p_src != NULL){ \
            __ml_node_##name *n = __ml_node_init_##name(); \
            if (n == NULL) return(-1); \
            n->data = p_src->data; \
            if (p_prev) p_prev->next = n; \
            else dest->head = n; \
            p_prev = n; \
            p_src = p_src->next; \
            ++n_e; \
        } \
        if (n_e != src->size){ \
            fprintf(stderr, "n_e=%u src->size=%u\n", n_e, src->size); \
            assert(n_e == src->size); \
        } \
        return(0); \
    } \
    \
    SCOPE int ml_insert_##name(ml_##name##_t *ll, ntype data, int skip_dup, int dup_ok){ \
        if (ll == NULL){ \
            fprintf(stderr, "error: ml_insert_%s: argument list is null\n", #name); \
            return(-1); \
        } \
        __ml_node_##name *n = __ml_node_init_##name(); \
        if (n == NULL) return(-1); \
        n->data = data; \
        \
        int cmp; \
        if (ll->head == NULL || ( (cmp = n_cmp(data, ll->head->data)) < 0 )){ \
            n->next = ll->head; \
            ll->head = n; \
            ++ll->size; \
            return(2); \
        } \
        \
        __ml_node_##name *p = ll->head; \
        while (cmp > 0 && p->next != NULL && ( (cmp = n_cmp(data, p->next->data)) > 0 )){ \
            p = p->next; \
        } \
        \
        if (cmp == 0 && dup_ok == 0){ \
            free(n); \
            fprintf(stderr, "error: ml_insert_%s: trying to add duplicate data\n", #name); \
            return(-1); \
        } \
        if (cmp == 0 && skip_dup != 0){ \
            free(n); \
            return(0); \
        } \
        \
        if (cmp <= 0 && p == ll->head){ \
            n->next = ll->head; \
            ll->head = n; \
        } else {\
            n->next = p->next; \
            p->next = n; \
        } \
        ++ll->size; \
        \
        if (cmp == 0) return(1); \
        else return(2); \
    } \
    \
    SCOPE __ml_node_##name *ml_head_##name(ml_##name##_t *ll){ \
        if (ll == NULL){ \
            fprintf(stderr, "error: ml_head_%s: argument list is null\n", #name); \
            return(NULL); \
        } \
        return(ll->head); \
    } \


/* types */
#define ml_t(name) ml_##name##_t
#define ml_node_t(name) __ml_node_##name

/* Declare data types and functions for generic linked list.
 * @param name type name of the linked list type
 * @param ntype type of the data value in node of the list.
 * @param cmp function that compares data values, returns -1 if x1 < x2, 0 if x1 == x2, 1 if x1 > x2
 * @return void
 */
#define ml_declare(name, ntype, cmp) ml_declare_i(static inline, name, ntype, cmp)

/* Allocate and initialize an empty ml list
 * @param ll pointer to propertly initialized list object from ml_alloc(name)
 * @return pointer to allocated list, or null on error.
 */
#define ml_alloc(name) ml_alloc_##name()
#define ml_init(name, ll) ml_init_##name(ll)

/* free or destroy a list
 * free does not free the list object @p ll itself, while dstry does.
 * @param name name of the macro list when declared, e.g. ml_declare(name, ...)
 * @param ll propertly initialized linked list object
 * @return void
 */
#define ml_free(name, ll) ml_free_##name(ll)
#define ml_dstry(name, ll) ml_dstry_##name(ll)

/* copy contents from src list to dest list.
 * @param src pointer to list of type ml_t(name)
 * @param dest pointer to list of type ml_t(name)
 * @return 0 on success, -1 on error
 */
#define ml_cpy(name, dest, src) ml_cpy_##name(dest, src)

/* ml_insert
 * @param name name of the macro list when declared
 * @param ll propertly initialized list object from ml_alloc(name)
 * @param data the data of type node to add, node is from declaration, e.g. ll_declare(name, node_type, ...)
 * @param skip_dup if non-zero and found, do nothing and return 0.
 * @param dup_ok if zero and duplicate found, print an error and return -1.
 * @return -1 on error, 0 if nothing added, 1 if duplicate added, 2 if new added
 */
#define ml_insert(name, ll, data, skip_dup, dup_ok) ml_insert_##name(ll, data, skip_dup, dup_ok)

// TODO: fix issue between invalid list and empty head
/* ml_head get head node of linked list
 * @param name name of the macro list when declared
 * @param ll pointer to propertly initialized list object from ml_alloc(name)
 * @return pointer to node of type ml_node_t(name)
 */
#define ml_begin(ll) ((ll)->head)

/* ml_node_next get next node in list
 * @param name name of the macro list when declared
 * @param node pointer to data of type ml_node_t(name)
 * @return next node in list, of type pointer to ml_node_t(name)
 */
#define ml_node_next(node) ((node)->next)

/* ml_node_val return the data value in a node
 * @param name name of the macro list when declared
 * @param node data of type ml_node_t(name)
 * @return data value of node
 */
#define ml_node_val(node) ((node)->data)

/* get number of elements in list
 * @param ll pointer to propertly initialized list object from ml_alloc(name)
 * @return uint32_t data
 */
#define ml_size(ll) ((ll)->size) 


/*******************************************************************************
 * vector
 ******************************************************************************/

#define mv_declare_i(SCOPE, name, type) \
    /* struct */ \
    typedef struct __mv_##name { \
        type *a; \
        uint32_t n, m; \
    } __mv_##name; \
    typedef struct __mv_##name mv_##name##_t; \
    \
    SCOPE int __mv_push_##name(mv_##name##_t *v, type x){ \
        if (v->n >= v->m){ \
            v->m = v->n + 1; \
            v->a = (type *)realloc(v->a, sizeof(type) * v->m); \
            if (v->a == NULL){ \
                fprintf(stderr, "error: __mv_push_%s: %s\n", #name, strerror(errno)); \
                return(-1); \
            } \
        } \
        \
        v->a[v->n++] = x; \
        return(0); \
    } \
    \
    SCOPE int __mv_insert_##name(mv_##name##_t *v, type x, uint32_t ix){ \
        if (ix >= v->m){ \
            v->m = ix + 1; \
            v->a = (type *)realloc(v->a, sizeof(type) * v->m); \
            if (v->a == NULL){ \
                fprintf(stderr, "error: __mv_insert_%s: %s\n", #name, strerror(errno)); \
                return(-1); \
            } \
        } \
        \
        v->a[ix] = x; \
        ++v->n; \
        return(0); \
    } \
    \
    SCOPE int __mv_resize_##name(mv_##name##_t *v, uint32_t n){ \
        if (n == v->m) return(0); \
        v->m = n; \
        v->a = (type *)realloc(v->a, sizeof(type) * v->m); \
        if (v->a == NULL){ \
            fprintf(stderr, "error: __mv_resize_%s: %s\n", #name, strerror(errno)); \
            return(-1); \
        } \
        return(0); \
    } \
    \

/* types */
#define mv_t(name) mv_##name##_t

/* declare */
#define mv_declare(name, type) mv_declare_i(static, name, type)

/* functions */

#define mv_init(v) ( (v)->n = (v)->m = 0, (v)->a = NULL )

#define mv_free(v) ( free((v)->a), (v)->a = NULL , (v)->n = (v)->m = 0 )

#define mv_i(v, ix) ((v)->a[(ix)])

#define mv_size(v) ( (v)->n )

#define mv_max(v) ( (v)->m )

#define mv_push(name, v, x) __mv_push_##name(v, x)

#define mv_insert(name, v, x, ix) __mv_insert_##name(v, x, ix)

#define mv_resize(name, v, n) __mv_resize_##name(v, n)

mv_declare(i32, int32_t);
mv_declare(u32, uint32_t);

/*******************************************************************************
 * binary tree
 ******************************************************************************/

/* key_cmp_fn is a function of type `int cmp(key_type key1, key_type key2)` that returns
* -1 if key1 < key2, 0 if key1 == key2, 1 if key1 > key2.
*  Values should be passed, not pointers to `key_type`
*/

#define mt_declare_i(SCOPE, name, key_type, value_type, key_cmp_fn) \
/* node struct*/ \
typedef struct __mt_node_##name { \
    key_type key; \
    value_type value; \
    KAVL_HEAD(struct __mt_node_##name) head; \
} __mt_node_##name; \
typedef struct __mt_node_##name mt_node_##name##_t; \
\
/* tree struct */ \
typedef struct __mt_tree_##name { \
    mt_node_##name##_t *root; \
} __mt_tree_##name; \
typedef struct __mt_tree_##name mt_tree_##name##_t; \
\
/* cmp function */ \
SCOPE int __key_cmp_##name(const mt_node_##name##_t *n1, const mt_node_##name##_t *n2){ \
    return key_cmp_fn(((n1)->key), ((n2)->key)); \
} \
\
KAVL_INIT2(__bt_##name, SCOPE, struct __mt_node_##name, head, __key_cmp_##name) \
\
/* tree iterator */ \
typedef struct __mt_itr_##name { \
    kavl_itr_t(__bt_##name) itr; \
    uint8_t has_data; \
} __mt_itr_##name; \
typedef struct __mt_itr_##name mt_itr_##name##_t; \
\
/* initialize tree */ \
SCOPE int __mt_tree_init_##name(__mt_tree_##name *b){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_tree_init_%s: tree is null\n", #name); \
        return(-1); \
    } \
    b->root = NULL; \
    return(0); \
} \
\
/* free tree */ \
SCOPE void __mt_tree_free_##name(__mt_tree_##name *b){ \
    if (b == NULL) \
        return; \
\
    if (b->root == NULL) \
        return; \
    kavl_itr_t(__bt_##name) itr; \
    kavl_itr_first(__bt_##name, b->root, &itr);  /* place at first */ \
    if (kavl_at(&itr) == NULL) \
        return; \
    do {                             /* traverse */ \
        mt_node_##name##_t *p = (mt_node_##name##_t *)kavl_at(&itr); \
        free((void*)p); \
    } while (kavl_itr_next(__bt_##name, &itr)); \
    b->root = NULL; \
\
}\
\
/* add to tree */\
SCOPE mt_node_##name##_t *__mt_tree_add_##name(__mt_tree_##name *b, \
    key_type key, value_type val, int *found){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_tree_add_%s: tree is null\n", #name); \
        return(NULL); \
    } \
\
    /* allocate node */ \
    mt_node_##name##_t *p = (mt_node_##name##_t *)calloc(1, sizeof(mt_node_##name##_t)); \
    if (p == NULL){ \
        fprintf(stderr, "error: __mt_tree_add_%s: %s\n", #name, strerror(errno)); \
        return(NULL); \
    } \
\
    /* initialize node */ \
    p->key = key; \
    p->value = val; \
\
    /* add node */ \
    mt_node_##name##_t *q; \
    q = kavl_insert(__bt_##name, &b->root, p, 0); \
    if (q != p) { \
        *found = 1; \
        free(p); \
    } else { \
        *found = 0; \
    } \
    return(q); \
}\
/* get node by key */ \
SCOPE mt_node_##name##_t *__mt_tree_find_##name(mt_tree_##name##_t *b, key_type key){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_tree_find_%s: tree is null\n", #name); \
        return(NULL); \
    } \
    mt_node_##name##_t p, *q; \
    p.key = key; \
    q = kavl_find(__bt_##name, b->root, &p, 0); \
    return q; \
} \
\
/* remove node */ \
SCOPE mt_node_##name##_t *__mt_tree_rm_##name(mt_tree_##name##_t *b, key_type key){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_tree_rm_%s: tree is null\n", #name); \
        return(NULL); \
    } \
    if (b->root == NULL) \
        return NULL; \
    mt_node_##name##_t p, *q; \
    p.key = key; \
    q = kavl_erase(__bt_##name, &b->root, &p, 0); \
    return q; \
} \
/* initialize iterator to first */ \
SCOPE int __mt_itr_first_##name(mt_tree_##name##_t *b, mt_itr_##name##_t *itr){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_itr_first_%s: tree is null\n", #name); \
        return(-1); \
    } \
    if (itr == NULL) {\
        fprintf(stderr, "error: __mt_itr_first_%s: iterator is null\n", #name); \
        return(-1); \
    } \
\
    if (b->root == NULL) {\
        itr->has_data = 0; \
        return 0; \
    } \
    kavl_itr_first(__bt_##name, b->root, &itr->itr); \
    if (kavl_at(&itr->itr) == NULL) { \
        itr->has_data = 0; \
    } else { \
        itr->has_data = 1; \
    } \
    return(0); \
} \
\
/* get next node in tree */ \
SCOPE int __mt_itr_next_##name(mt_itr_##name##_t *itr){ \
    if (itr == NULL) {\
        fprintf(stderr, "error: __mt_itr_next_%s: iterator is null\n", #name); \
        return(-1); \
    } \
    if (itr->has_data == 0) \
        return 0; \
    /* calling kavl_itr_next places the iterator to point to the next node in the tree \
     * if this next node is valid and contains data, 1 is returned and kavl_at(itr) returns the data. \
     * if this next node is empty, 0 is returned and kavl_at(itr) returns NULL. \
     */ \
    itr->has_data = kavl_itr_next(__bt_##name, &itr->itr); \
    return(itr->has_data); \
} \
\
/* check whether iterator contains data */ \
SCOPE int __mt_itr_valid_##name(mt_itr_##name##_t *itr){ \
    if (itr == NULL) \
        return 0; \
    return itr->has_data; \
} \
\
/* place iterator at key value */ \
SCOPE int __mt_itr_find_##name(mt_tree_##name##_t *b, key_type key, mt_itr_##name##_t *itr){ \
    if (b == NULL) {\
        fprintf(stderr, "error: __mt_itr_find_%s: tree is null\n", #name); \
        return(-1); \
    } \
    if (itr == NULL) {\
        fprintf(stderr, "error: __mt_itr_find_%s: iterator is null\n", #name); \
        return(-1); \
    } \
    if (b->root == NULL) {\
        return 0; \
    } \
\
    mt_node_##name##_t p; \
    p.key = key; \
    int ret = kavl_itr_find(__bt_##name, b->root, &p, &itr->itr); \
\
    if (kavl_at(&itr->itr) == NULL) { \
        itr->has_data = 0; \
    } else { \
        itr->has_data = 1; \
    } \
\
    return(ret); \
} \
\
/* iterator key */ \
SCOPE key_type *__mt_itr_key_##name(mt_itr_##name##_t *itr){ \
    if (itr == NULL) { \
        fprintf(stderr, "error: __mt_itr_key_%s: iterator is null\n", #name); \
        return(NULL); \
    } \
    mt_node_##name##_t *p = (mt_node_##name##_t *)kavl_at(&itr->itr); \
    if (p == NULL) \
        return NULL; \
    return(&p->key); \
} \
\
/* iterator value */ \
SCOPE value_type *__mt_itr_val_##name(mt_itr_##name##_t *itr){ \
    if (itr == NULL) {\
        fprintf(stderr, "error: __mt_itr_val_%s: iterator is null\n", #name); \
        return(NULL); \
    } \
    mt_node_##name##_t *p = (mt_node_##name##_t *)kavl_at(&itr->itr); \
    if (p == NULL) \
        return NULL; \
    return &p->value; \
}

/* declare type and functions with static scope */
#define mt_declare(name, key_type, value_type, key_cmp_fn) \
    mt_declare_i(static inline, name, key_type, value_type, key_cmp_fn)

/* types
 * @param name The name specifier of the tree type
 */
#define mt_node_t(name) mt_node_##name##_t
#define mt_tree_t(name) mt_tree_##name##_t
#define mt_itr_t(name) mt_itr_##name##_t

/* functions */

/* initialize tree
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 *
 * @retval `int`
 * @return 0 on success, -1 on error
 */
#define mt_init(name, tree) __mt_tree_init_##name((tree))

/* free nodes in tree
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 *
 * @return void
 */
#define mt_free(name, tree) __mt_tree_free_##name((tree))

/* add node to tree
 *
 * Dynamically allocates a node with key_type and value_type, and assigns
 * `key` and `val`by simple copy assignment.
 *
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 * @param key key of type `key_type`
 * @param val value of type `value_type`
 * @param found pointer to int to hold return value of found status
 *
 * @retval `mt_node_t(name) *`
 * @return Pointer to mt_node_t(name) with key-value data.
 *         If key already exists, returns pointer to existing node and sets *found=1.
 *         Otherwise, returns pointer to newly allocated node and sets *found=0.
 *         Will return NULL on error.
 * @note Unclear exactly how kavl_insert will return on error.
 */
#define mt_add(name, tree, key, val, found) __mt_tree_add_##name((tree), (key), (val), (found))

/* find node by key value
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 * @param key key of type `key_type`
 *
 * @retval `mt_node_t(name) *`
 * @return Pointer to node of type mt_node_t(name) with key-value data.
 *         If key is not found, returns NULL.
 */
#define mt_find(name, tree, key) __mt_tree_find_##name((tree), (key))

/* remove node by key value
 * After last element is removed, `tree.root` is set to NULL.
 *
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 * @param key key of type `key_type`
 *
 * @retval `mt_node_t(name) *`
 * @return Pointer to removed node of type mt_node_t(name) with key-value data.
 *         If key is not found, returns NULL.
 */
#define mt_remove(name, tree, key) __mt_tree_rm_##name((tree), (key))

/* initialize iterator to first node in tree
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `int`
 * @return 0 on success, -1 on error
 */
#define mt_itr_first(name, tree, itr) __mt_itr_first_##name((tree), (itr))

/* get next node in tree
 * @param name tree type specifier
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `int`
 * @return 0 if no next element, 1 if next element found
 */
#define mt_itr_next(name, itr) __mt_itr_next_##name((itr))

/* find node by key value
 *
 * places iterator at node with key value. if key is not found, iterator
 *  1) is NULL if key is greater than all keys in the tree.
 *  2) points to the next node greater than key if it exists.
 *
 * @param name tree type specifier
 * @param tree pointer to mt_tree_t(name)
 * @param key key of type `key_type`
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `int`
 * @return 0 if key not found, 1 if key found
 */
#define mt_itr_find(name, tree, key, itr) __mt_itr_find_##name((tree), (key), (itr))

/* check if iterator has data.
 * @param name tree type specifier
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `int`
 * @return 0 if iterator has no data, 1 if iterator has data
 */
#define mt_itr_valid(name, itr) __mt_itr_valid_##name((itr))

/* return the key of the node at the iterator
 * @param name tree type specifier
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `value_type *`
 * @return Pointer to key in tree. NULL if itr has no data.
 */
#define mt_itr_key(name, itr) __mt_itr_key_##name((itr))

/* return the value of the node at the iterator
 * @param name tree type specifier
 * @param itr pointer to mt_itr_t(name)
 *
 * @retval `value_type *`
 * @return Pointer to value in tree. NULL if itr has no data.
 */
#define mt_itr_val(name, itr) __mt_itr_val_##name((itr))

/* return number of elements in tree
 * @param tree pointer to mt_tree_t(name)
 *
 * @retval `uint32_t`
 * @return number of elements
 */
#define mt_size(tree) kavl_size(head, ((tree)->root))

/* Example
#include "g_list.h"
#include <stdio.h>
#include <string.h>

// Define key and value types
typedef char *key_type;
typedef int value_type;

// Comparison function for keys
int key_cmp(key_type k1, key_type k2) { return strcmp(k1, k2); }

// Declare the mt_tree
mt_declare_i(static, str_int, key_type, value_type, key_cmp)

int main() {
    // Initialize the tree
    mt_tree_t(str_int) tree;
    mt_init(str_int, &tree);

    // Add some key-value pairs
    char *keys[] = {"apple", "banana", "cherry", "apple", "peach"};
    int values[] = {1, 2, 3, 4, 5};
    mt_node_t(str_int) * node;
    for (int i = 0; i < 5; i++) {
        node = mt_add(str_int, &tree, keys[i], values[i]);
        printf("Added key '%s' with value %d\n", node->key, node->value);
    }

    // Iterate over the tree
    int iret;
    mt_itr_t(str_int) itr;
    mt_itr_first(str_int, &tree, &itr);
    printf("\nBag contents:\n");
    while (mt_itr_valid(str_int, &itr)) {
        key_type *pkey = mt_itr_key(str_int, &itr);
        value_type *value = mt_itr_val(str_int, &itr);
        printf("Key: '%s', Value: %d, next: %i\n", *pkey, *value, iret);
        iret = mt_itr_next(str_int, &itr);
    }

    // find with iter
    key_type key = "banana";
    printf("\nFinding key '%s' by iterator\n", key);
    int fret = mt_itr_find(str_int, &tree, key, &itr);
    if (fret == 1) {
        printf("Found key '%s' with value %d\n", *mt_itr_key(str_int, &itr),
               *mt_itr_val(str_int, &itr));
    } else {
        printf("Key '%s' not found\n", key);
    }

    // perform a find operation
    printf("\nFinding key '%s' by node\n", key);
    node = mt_find(str_int, &tree, key);
    if (node) {
        printf("Found key '%s' with value %d\n", node->key, node->value);
    } else {
        printf("Key '%s' not found\n", key);
    }

    // remove banana
    printf("\nRemoving key '%s'\n", key);
    key = "banana";
    node = mt_remove(str_int, &tree, key);
    if (node)
        printf("Removed key '%s' with value %d\n", node->key, node->value);
    node = mt_find(str_int, &tree, key);
    if (node) {
        printf("Key '%s' found with value %d\n", node->key, node->value);
    } else {
        printf("Key '%s' not found\n", key);
    }

    // Free the tree
    mt_free(str_int, &tree);

    return 0;
}
*/

#endif // G_LIST_H
