/*
 * Allocate and deallocate BTN data structures and BTN trees.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <stdlib.h>		/* For `malloc()' and `free()' */
#include "fntype.h"
#include "btn.h"

/***************************************************************************
 * makebtn()
 *
 *	Allocate a BTN data structure with given members.
 *
 *  Calling sequence:
 * 	makebtn( content, left, right, tag )
 *
 *  Basic algorithm:
 *
 *   Allocate a BTN data structure at NODE
 *   Let NODE.TAG = TAG
 *   Let NODE.CONTENT  = CONTENT
 *   Let NODE.LEFT = LEFT
 *   Let NODE.RIGHT = RIGHT
 *   Return NODE
 *
 *
 *  Inputs:
 *      (void *)content		This can be `&interval', `&real' or NULL.
 * 
 * 	(btn *)left		These can be valid BTN pointers, or NULL.
 *	(btn *)right
 * 
 *	(void *)tag		This points to a data structure, or is NULL.
 *
 *  Outputs:
 *	(btn *)makebtn		The return value is a pointer to a newly
 *				  allocated BTN structure.
 *
 *
 * External functions called:
 *	malloc()	Declared in <stdlib.h>
 */
extern btn *
  makebtn(
	  void *content,	/* Content data structure, or NULL. */
	  btn  *left,		/* Left descendent node, or NULL. */
	  btn  *right,		/* Right descendent node, or NULL. */
	  void *tag)		/* Extra information, or NULL. */
{
  btn *node;

  node = (btn *)malloc(sizeof(btn));  assert(node);
  node->content = content;
  node->left = left;
  node->right = right;
  node->tag = tag;

  return(node);
}

/***************************************************************************
 * freebtn()
 *
 *	Deallocate a BTN data structure and any non-NULL members.
 *
 *  Calling sequence and basic algorithm:
 *
 * 	freebtn( NODE, FREECONTENT, FREETAG ):
 *         If NODE != NULL then
 *            If NODE.CONTENT != NULL && FREECONTENT != NULL then
 *               Deallocate NODE.CONTENT with FREECONTENT()
 *            If NODE.TAG != NULL && FREETAG != NULL then
 *               Deallocate NODE.TAG with FREETAG
 *            Deallocate NODE
 *         Return NULL
 *
 *
 *  Input:
 * 	(btn *)node		This is a pre-allocated BTN data structure.
 *
 *	(freetype)freecontent	This function will be used to deallocate
 *				  the content member.
 *
 *	(freetype)freetag	This function will be used to deallocate
 *				  the tag member.
 *
 *  Outputs:
 *	The return value is always a NULL pointer.
 *
 *
 * External functions called:
 *	free()	Declared in <stdlib.h>
 */
extern btn *
  freebtn(
	  btn  *node,		/* Previously-allocated BTN structure. */
	  freetype freecontent, /* Used to free the content member.      */
	  freetype freetag)	/* Used to free the tag member.          */
{
  if(node)
    {
      if(node->content && freecontent)
	node->content = freecontent(node->content);
      if(node->tag && freetag)
	node->tag = freetag(node->tag);
      free(node);
    }
  return(0);
}

/***************************************************************************
 * makebtnt()
 *
 *	Allocate a complete binary tree of BTN data structures down
 *	to a given level.  All `content' and `tag' members are NULL.
 *
 *  Calling sequence:
 * 	makebtnt( level )
 *
 *  Basic algorithm:
 *
 *   Let ROOT = makebtn(NULL,NULL,NULL,NULL)
 *   If LEVEL>0 then
 *      Let ROOT.LEFT = makebtnt( LEVEL-1 )
 *      Let ROOT.RIGHT = makebtnt( LEVEL-1 )
 *   Return ROOT
 *
 *
 *  Inputs:
 *	(int)level	This should be a small nonnegative integer.
 *
 *  Outputs:
 *	(btn *)makebtnt	The return value is a pointer to the root of
 *			a newly allocated tree of BTN structures.
 *
 *  External function called:
 *	makebtn()
 *
 *  Assumptions:
 *	`level>=0' and `level<32'.
 */
extern btn *
  makebtnt(
	   int level)		/* Number of levels in the tree. */
{
  btn *root;

  assert(level>=0);
  assert(level<32);		/* Trees grow exponentially with depth. */

  root = makebtn(0,0,0,0); assert(root);
  if(level>0)
    {
      root->left  = makebtnt(level-1);
      root->right = makebtnt(level-1);
    }
  return(root);
}

/***************************************************************************
 * freebtnt()
 *
 *	Deallocate a tree of BTN data structures and any non-NULL members.
 *
 *  Calling sequence and basic algorithm:
 *
 * 	freebtnt( ROOT, FREECONTENT, FREETAG ):
 *         If ROOT != NULL then
 *            freebtnt( ROOT.LEFT, FREECONTENT, FREETAG )
 *            freebtnt( ROOT.RIGHT, FREECONTENT, FREETAG )
 *            freebtn( ROOT, FREECONTENT, FREETAG )
 *         Return NULL
 *      
 *
 *  Input:
 * 	(btn *)root	This should be a previously allocated and linked
 *			binary tree of BTN data structures.
 *
 *	(freetype)freecontent	This function will be used to deallocate
 *				  the content member.
 *
 *	(freetype)freetag	This function will be used to deallocate
 *				  the tag member.
 *
 *  Outputs:
 *	The return value is always a NULL pointer.
 *
 * External functions called:
 *	free()		Declared in <stdlib.h>
 *	freebtn()
 *
 *  Assumption:
 *	Leaf BTNs have NULL in both their `left' and `right' members.
 *
 */
extern btn *
  freebtnt(
	   btn  *root,		/* Root of pre-allocated BTN tree. */
	   freetype freecontent, /* Used to free the content member. */
	   freetype freetag)	/* Used to free the tag member.      */
{
  if(root)
    {
      freebtnt(root->left, freecontent, freetag);
      freebtnt(root->right, freecontent, freetag);
      freebtn(root, freecontent, freetag);
    }
  return(0);
}

/***************************************************************************
 * btn2branch()
 *
 *	Allocate a BTN tree branch out to a target node.
 *
 *  Calling sequence:
 * 	btn2branch( self, level, block )
 *
 *  Basic algorithm:
 *
 *  If LEVEL>0 then
 *    If BLOCK is even then
 *      If SELF.LEFT==NULL then
 *        Let SELF.LEFT = makebtn( NULL, NULL, NULL, NULL )
 *      Let SELF = btn2branch( SELF.LEFT, LEVEL-1, BLOCK/2 )
 *    Else
 *      If SELF.RIGHT==NULL then
 *        Let SELF.RIGHT = makebtn( NULL, NULL, NULL, NULL )
 *      Let SELF = btn2branch( SELF.RIGHT, LEVEL-1, (BLOCK-1)/2 )
 *  Return SELF
 *
 *
 *  Input:
 * 	(btn *)self	This should be a previously allocated BTN.
 *
 * 	(int)level	This should be a small nonnegative integer.
 *
 * 	(int)block	This should be a nonnegative integer.
 *
 *  Outputs:
 *	(btn *)btn2branch	The return value is the leaf at 
 *				 (`level',`block') in a partial
 *				 BTN tree.  The leaf is created
 *				 if necessary.
 *
 * External functions called:
 *	makebtn()
 *
 *  Assumptions:
 *	1.   `level==0' is the root of the tree.
 *	2.   `block' is interpreted in Paley or natural order:
 *		left branch is even, right branch is odd.
 *	3.   `self' is non-NULL.
 *	4.   `level<256' 
 *	5.   `block>=0'
 *
 */
extern btn *
  btn2branch(
	     btn  *self,	/* Current root of the BTN tree. */
	     int  level,	/* Target level with respect to `self'. */
	     int  block)	/* Block index with respect to `self'. */
{
  assert(level>=0);
  assert(level<256);
  assert(self);
  assert(block>=0);

  if(level>0)
    {
      if(block&1)			/* `block' is odd  ==> branch to right. */
	{
	  if( self->right == 0 )
	    self->right = makebtn(0,0,0,0);
	  self = btn2branch( self->right, level-1, (block-1)/2 );
	}
      else		/* `block' is even ==> branch to left. */
	{
	  if( self->left == 0 )
	    self->left = makebtn(0,0,0,0);
	  self = btn2branch( self->left, level-1, block/2 );
	}
    }
  return(self);
}

/***************************************************************************
 * btnt2btn()
 *
 *	Get a BTN from a binary tree of BTNs
 *
 *  Calling sequence:
 * 	btnt2btn( root, level, block )
 *
 *  Basic algorithm:
 *
 *   If LEVEL==0 || ROOT==NULL then
 *      Let NODE = ROOT
 *   Else
 *      If BLOCK is even then
 *         Let NODE = btnt2btn( ROOT.LEFT, LEVEL-1, BLOCK/2 )
 *      Else
 *         Let NODE = btnt2btn( ROOT.RIGHT, LEVEL-1, (BLOCK-1)/2 )
 *   Return NODE
 *
 *
 *  Input:
 * 	(btn *)root	This should be a previously allocated BTN.
 *
 * 	(int)level	This should be a small nonnegative integer.
 *
 * 	(int)block	This should be a nonnegative integer.
 *
 *  Outputs:
 *	(btn *)btnt2btn	The return value is the leaf at (`level',`block')
 *			 in a partial BTN tree.  If the tree ends before
 *			 the leaf, NULL is returned.
 *
 *  Assumptions:
 *	1.   `level==0' is the root of the tree.
 *	2.   `block' is interpreted in Paley or natural order:
 *		left branch is even, right branch is odd.
 *	3.   `root' is non-NULL.
 *	4.   `level<256' 
 *	5.   `block>=0'
 *
 */
extern btn *
  btnt2btn(
	   btn  *root,		/* Current root of the BTN tree. */
	   int  level,		/* Target level with respect to `root'. */
	   int  block)		/* Block index with respect to `root'. */
{
  btn *node;

  assert(level>=0);
  assert(level<256);
  assert(root);
  assert(block>=0);

  if(level==0 || root==0)
    node = root;
  else
    {
      if(block&1)		/* `block' is odd  ==> branch to right. */
	node = btn2branch( root->right, level-1, (block-1)/2 );
      else			/* `block' is even ==> branch to left. */
	node = btnt2btn( root->left, level-1, block/2 );
    }
  return(node);
}
