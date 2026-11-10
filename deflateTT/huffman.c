/*
 * Huffman codes for deflateTT.
 * Programmed by Ville Helin in 2007, with length limiting added in 2026.
 *
 * This code is under GNU General Public Licence (GPL), version 2, June 1991.
 *
 * A symbol that is alone in its alphabet would otherwise get length 0,
 * which the decoder cannot read. The same is true of an empty alphabet.
 * Both cases are given two length-1 codes, the way zlib does for the
 * distance tree. Lengths longer than 15 are folded back with the same
 * adjustment zlib uses in gen_bitlen, then canonical codes are assigned
 * as in RFC 1951 section 3.2.2.
 */

#include <stdlib.h>

#include "huffman.h"

#define HUFFMAN_MAX_BITS 15

struct node {
  struct node *left;
  struct node *right;
  struct node *smaller;
  struct node *larger;
  int literal;
  int weight;
  int leaf;
  int code;
  int codeLength;
};

static void pq_push(struct node **queue, struct node *item) {

  struct node *q;

  item->smaller = NULL;
  item->larger = NULL;
  q = *queue;

  if (q == NULL) {
    *queue = item;
    return;
  }

  while (1) {
    if (item->weight > q->weight) {
      if (q->larger != NULL)
        q = q->larger;
      else {
        q->larger = item;
        item->smaller = q;
        return;
      }
    }
    else {
      item->smaller = q->smaller;
      item->larger = q;
      if (q->smaller != NULL)
        q->smaller->larger = item;
      q->smaller = item;

      if (q == *queue)
        *queue = item;

      return;
    }
  }
}


static struct node *pq_pop(struct node **queue) {

  struct node *item;

  item = *queue;
  if (item == NULL)
    return NULL;

  *queue = item->larger;
  if (*queue != NULL)
    (*queue)->smaller = NULL;

  return item;
}


static void propagate_codes(struct node *item, int code, int codeLength) {

  if (item->leaf) {
    item->code = code;
    item->codeLength = codeLength;
    return;
  }

  if (item->left != NULL)
    propagate_codes(item->left, (code << 1) | 1, codeLength + 1);
  if (item->right != NULL)
    propagate_codes(item->right, (code << 1) | 0, codeLength + 1);
}


static void free_parents(struct node **parents, int parentCount) {

  int i;

  for (i = 0; i < parentCount; i++)
    free(parents[i]);
}


/* Give every alphabet at least two positive weights so the tree has a
   real bit and never a length-0 code. Dummy symbols are not emitted. */
static int ensure_two_symbols(struct node *nodes, int n) {

  int i, used, dummy;

  used = 0;
  dummy = -1;

  for (i = 0; i < n; i++) {
    if (nodes[i].weight > 0)
      used++;
    else if (dummy < 0)
      dummy = i;
  }

  if (used >= 2)
    return 0;

  if (used == 0) {
    if (n < 2)
      return -1;
    nodes[0].weight = 1;
    nodes[1].weight = 1;
    return 0;
  }

  if (dummy < 0)
    return -1;

  nodes[dummy].weight = 1;
  return 0;
}


static void sort_by_freq(int *syms, int count, const int *freq) {

  int a, b, key;

  for (a = 1; a < count; a++) {
    key = syms[a];
    b = a;
    while (b > 0 && (freq[syms[b - 1]] > freq[key] || (freq[syms[b - 1]] == freq[key] && syms[b - 1] > key))) {
      syms[b] = syms[b - 1];
      b--;
    }
    syms[b] = key;
  }
}


/* Count leaf lengths, capping the walk the way zlib's gen_bitlen does.
   Internal nodes past maxLen add to overflow too. That extra count is
   what makes the following adjustment restore a complete code. */
static void accumulate_lengths(struct node *node, int depth, int maxLen,
                               int *count, int *overflow) {

  int childDepth;

  if (node->leaf) {
    if (depth > maxLen) {
      (*overflow)++;
      count[maxLen]++;
    }
    else if (depth > 0)
      count[depth]++;
    return;
  }

  if (depth > maxLen)
    (*overflow)++;

  childDepth = (depth > maxLen ? maxLen : depth) + 1;
  if (node->left != NULL)
    accumulate_lengths(node->left, childDepth, maxLen, count, overflow);
  if (node->right != NULL)
    accumulate_lengths(node->right, childDepth, maxLen, count, overflow);
}


static int limit_from_tree(struct node *tree, struct node *nodes,
                           const int *freq, int n, int *lengths) {

  int count[HUFFMAN_MAX_BITS + 1];
  int i, bits, overflow, used, idx, left;
  int *syms;

  overflow = 0;
  for (i = 0; i <= HUFFMAN_MAX_BITS; i++)
    count[i] = 0;

  accumulate_lengths(tree, 0, HUFFMAN_MAX_BITS, count, &overflow);

  if (overflow == 0) {
    for (i = 0; i < n; i++) {
      lengths[i] = nodes[i].codeLength;
      if (lengths[i] > HUFFMAN_MAX_BITS)
        return -1;
    }
    return 0;
  }

  used = 0;
  for (i = 0; i < n; i++) {
    if (nodes[i].weight > 0)
      used++;
  }
  if (used > (1 << HUFFMAN_MAX_BITS))
    return -1;

  while (overflow > 0) {
    bits = HUFFMAN_MAX_BITS - 1;
    while (bits > 0 && count[bits] == 0)
      bits--;
    if (bits <= 0)
      return -1;
    count[bits]--;
    count[bits + 1] += 2;
    if (count[HUFFMAN_MAX_BITS] <= 0)
      return -1;
    count[HUFFMAN_MAX_BITS]--;
    overflow -= 2;
  }

  syms = (int *)malloc(sizeof(int) * used);
  if (syms == NULL)
    return -1;

  idx = 0;
  for (i = 0; i < n; i++) {
    lengths[i] = 0;
    if (nodes[i].weight > 0)
      syms[idx++] = i;
  }
  if (idx != used) {
    free(syms);
    return -1;
  }

  sort_by_freq(syms, used, freq);

  idx = 0;
  for (bits = HUFFMAN_MAX_BITS; bits >= 1; bits--) {
    int leftAtLength;

    leftAtLength = count[bits];
    while (leftAtLength > 0) {
      if (idx >= used) {
        free(syms);
        return -1;
      }
      lengths[syms[idx]] = bits;
      idx++;
      leftAtLength--;
    }
  }
  free(syms);

  if (idx != used)
    return -1;

  left = 1;
  for (bits = 1; bits <= HUFFMAN_MAX_BITS; bits++) {
    left <<= 1;
    left -= count[bits];
    if (left < 0)
      return -1;
  }

  return 0;
}


static int assign_canonical(int n, const int *lengths, int *codes) {

  int count[HUFFMAN_MAX_BITS + 1];
  int next[HUFFMAN_MAX_BITS + 1];
  int i, bits, length, code, left;

  for (i = 0; i <= HUFFMAN_MAX_BITS; i++)
    count[i] = 0;

  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length < 0 || length > HUFFMAN_MAX_BITS)
      return -1;
    count[length]++;
  }

  left = 1;
  for (bits = 1; bits <= HUFFMAN_MAX_BITS; bits++) {
    left <<= 1;
    left -= count[bits];
    if (left < 0)
      return -1;
  }

  code = 0;
  count[0] = 0;
  for (bits = 1; bits <= HUFFMAN_MAX_BITS; bits++) {
    code = (code + count[bits - 1]) << 1;
    next[bits] = code;
  }

  for (i = 0; i < n; i++) {
    length = lengths[i];
    if (length != 0) {
      codes[i] = next[length];
      next[length]++;
    }
    else
      codes[i] = 0;
  }

  return 0;
}


int huffmanBuildCodes(const int *frequencies, int n, int *codeLengths, int *codes) {

  struct node *nodes, *tree, *node1, *node2, *parent;
  struct node **parents, *queue;
  int i, parentCount;

  if (frequencies == NULL || codeLengths == NULL || codes == NULL)
    return -1;
  if (n < 1 || n > 4096)
    return -1;

  for (i = 0; i < n; i++) {
    if (frequencies[i] < 0)
      return -1;
  }

  if (n == 1) {
    codeLengths[0] = frequencies[0] > 0 ? 1 : 0;
    return assign_canonical(1, codeLengths, codes);
  }

  nodes = (struct node *)malloc(sizeof(struct node) * n);
  parents = (struct node **)malloc(sizeof(struct node *) * n);
  if (nodes == NULL || parents == NULL) {
    free(nodes);
    free(parents);
    return -1;
  }

  parentCount = 0;
  queue = NULL;
  tree = NULL;

  for (i = 0; i < n; i++) {
    nodes[i].left = NULL;
    nodes[i].right = NULL;
    nodes[i].smaller = NULL;
    nodes[i].larger = NULL;
    nodes[i].literal = i;
    nodes[i].weight = frequencies[i];
    nodes[i].leaf = 1;
    nodes[i].code = 0;
    nodes[i].codeLength = 0;
  }

  if (ensure_two_symbols(nodes, n) != 0) {
    free_parents(parents, parentCount);
    free(parents);
    free(nodes);
    return -1;
  }

  for (i = 0; i < n; i++) {
    if (nodes[i].weight > 0)
      pq_push(&queue, &nodes[i]);
  }

  while (1) {
    node1 = pq_pop(&queue);
    node2 = pq_pop(&queue);

    if (node2 == NULL) {
      tree = node1;
      break;
    }

    parent = (struct node *)malloc(sizeof(struct node));
    if (parent == NULL || parentCount >= n) {
      free(parent);
      free_parents(parents, parentCount);
      free(parents);
      free(nodes);
      return -1;
    }

    parent->code = 0;
    parent->codeLength = 0;
    parent->leaf = 0;
    parent->literal = -1;
    parent->weight = node1->weight + node2->weight;
    parent->left = node1;
    parent->right = node2;
    parent->smaller = NULL;
    parent->larger = NULL;
    parents[parentCount++] = parent;
    pq_push(&queue, parent);
  }

  if (tree == NULL) {
    free_parents(parents, parentCount);
    free(parents);
    free(nodes);
    return -1;
  }

  propagate_codes(tree, 0, 0);

  if (limit_from_tree(tree, nodes, frequencies, n, codeLengths) != 0) {
    free_parents(parents, parentCount);
    free(parents);
    free(nodes);
    return -1;
  }

  free_parents(parents, parentCount);
  free(parents);
  free(nodes);

  return assign_canonical(n, codeLengths, codes);
}
