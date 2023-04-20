
#ifndef G_LIST_H
#define G_LIST_H

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>

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
                fprintf(stderr, "error: __mv_push_%s: %s\n", #name, strerror(errno)); \
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

#define mv_i(v, ix) (v)->a[(ix)]

#define mv_size(v) ( (v)->n )

#define mv_max(v) ( (v)->m )

#define mv_push(name, v, x) __mv_push_##name(v, x)

#define mv_insert(name, v, x, ix) __mv_insert_##name(v, x, ix)

#define mv_resize(name, v, n) __mv_resize_##name(v, n)

mv_declare(i32, int32_t);

#endif // G_LIST_H
