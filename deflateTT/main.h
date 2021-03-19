
#ifndef _MAIN_H
#define _MAIN_H

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

#endif
