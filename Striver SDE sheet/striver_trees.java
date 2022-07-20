//! Striver SDE sheet all Trees solution

import java.util.*;

public class striver_trees {
    // ~ Binary Tree definition
    public class TreeNode {
        int val;
        TreeNode left;
        TreeNode right;

        TreeNode() {
        }

        TreeNode(int val) {
            this.val = val;
        }

        TreeNode(int val, TreeNode left, TreeNode right) {
            this.val = val;
            this.left = left;
            this.right = right;
        }
    }

    // ? Day 17 Binary Tree

    // ^ Binary Tree Inorder Traversal
    // Given the root of a binary tree, return the inorder traversal of its nodes'
    // values.

    class inorder_recursion {

        // declare the list outside the function, else it won't work
        List<Integer> list = new ArrayList<>();

        public List<Integer> inorderTraversal(TreeNode root) {
            if (root == null)
                return list;

            inorderTraversal(root.left);
            list.add(root.val);
            inorderTraversal(root.right);

            return list;
        }
    }

    // ^ Binary Tree Preorder Traversal

    // Given the root of a binary tree, return the preorder traversal of its nodes'
    // values.

    class preorder_recursive {

        // declare resultant list as global variable
        List<Integer> list = new ArrayList<>();

        public List<Integer> preorderTraversal(TreeNode root) {
            if (root == null)
                return list;

            list.add(root.val);
            preorderTraversal(root.left);
            preorderTraversal(root.right);

            return list;
        }
    }

    // ^Binary Tree Postorder Traversal

    // Given the root of a binary tree, return the postorder traversal of its nodes'
    // values.

    class postorder_traversal {

        // declare resultant array list as global variable
        List<Integer> list = new ArrayList<>();

        public List<Integer> postorderTraversal(TreeNode root) {
            if (root == null)
                return list;

            postorderTraversal(root.left);
            postorderTraversal(root.right);
            list.add(root.val);

            return list;
        }
    }

    // ^ Morris inorder Traversal

    // Given the root of a binary tree, return the inorder traversal of its nodes'
    // values using iterative method

    class inorder_iterative {
        public List<Integer> inorderTraversal(TreeNode root) {
            List<Integer> list = new ArrayList<>();
            Stack<TreeNode> s = new Stack<>();
            TreeNode node = root; // a node to traverse

            // go to left node first
            while (true) {
                // if node is not null, we push it to the stack and move to the left
                if (node != null) {
                    s.push(node);
                    node = node.left;
                } else { // when it is null
                    if (s.isEmpty()) {
                        // break when the stack goes empty
                        break;
                    }
                    // or add the top element of the stack and move to the right
                    node = s.pop();
                    list.add(node.val);
                    node = node.right;
                }
            }
            return list;
        }
    }

    // ^ Morris Preorder traversal

    // ^ Left View of Binary Tree

    // ^ Bottom View of Binary Tree

    // ^ Top View of Binary Tree

    // ^ Preorder inorder postorder in a single traversal

    // ^ Vertical order traversal

    // ^ Root to node path in a Binary Tree

    // ^ Max width of a Binary Tree

    // ^ Level order Traversal / Level order traversal in spiral form

    // Given the root of a binary tree, return the level order traversal of its
    // nodes' values. (i.e., from left to right, level by level).

    // Input: root = [3,9,20,null,null,15,7]
    // Output: [[3],[9,20],[15,7]]

    class level_order {
        public List<List<Integer>> levelOrder(TreeNode root) {
            Queue<TreeNode> q = new LinkedList<>();
            List<List<Integer>> ans = new LinkedList<List<Integer>>();

            if (root == null)
                return ans;
            // first insert the root in the queue
            q.offer(root);
            while (!q.isEmpty()) {
                int levelSize = q.size();// stores current queue size
                List<Integer> list = new LinkedList<>();
                for (int i = 0; i < levelSize; i++) { // iterate for the current queue size
                    // if left is not null of current queue element, we add it to the queue
                    if (q.peek().left != null)
                        q.offer(q.peek().left);
                    // similarly for right
                    if (q.peek().right != null)
                        q.offer(q.peek().right);

                    // add the current queue value to the list
                    list.add(q.poll().val);
                }
                // add the list to the resultant list
                ans.add(list);
            }
            return ans;
        }
    }

    // ^ Height of a Binary Tree

    // Given the root of a binary tree, return its maximum depth.

    // A binary tree's maximum depth is the number of nodes along the longest path
    // from the root node down to the farthest leaf node.

    // Example 1:

    // Input: root = [3,9,20,null,null,15,7]
    // Output: 3

    class tree_height {
        public int maxDepth(TreeNode root) {
            if (root == null)
                return 0;

            int left = maxDepth(root.left);
            int right = maxDepth(root.right);

            return Math.max(left, right) + 1;
        }
    }

    // ^ Diameter of Binary Tree

    // Given the root of a binary tree, return the length of the diameter of the
    // tree.

    // The diameter of a binary tree is the length of the longest path between any
    // two nodes in a tree. This path may or may not pass through the root.

    // The length of a path between two nodes is represented by the number of edges
    // between them.

    // Example 1:

    // Input: root = [1,2,3,4,5]
    // Output: 3
    // Explanation: 3 is the length of the path [4,2,1,3] or [5,2,1,3].

    class diameter_of_tree {
        public int diameterOfBinaryTree(TreeNode root) {
            // just create a 1 size array to store the answer in a smart way
            int[] diameter = new int[1];
            height(root, diameter);
            return diameter[0];
        }

        // just normal height function with 1 extra line to find max of left height and
        // right height
        private int height(TreeNode root, int[] diameter) {
            if (root == null)
                return 0;

            int lh = height(root.left, diameter);
            int rh = height(root.right, diameter);

            // just store the maximum sum of left and right subtree, thats our answer
            diameter[0] = Math.max(lh + rh, diameter[0]);

            return 1 + Math.max(lh, rh);
        }
    }

    // ^ Check if the Binary tree is height-balanced or not

    // Given a binary tree, determine if it is height-balanced.

    // For this problem, a height-balanced binary tree is defined as:

    // a binary tree in which the left and right subtrees of every node differ in
    // height by no more than 1.

    class height_balanced {
        public boolean isBalanced(TreeNode root) {
            return height(root) != -1;
        }

        // simple height calculation method with a few changes
        public int height(TreeNode root) {
            if (root == null)
                return 0;

            // finding height of left subtree and right subtree
            int leftHeight = height(root.left);
            // if at any point, the left or right subtree height is -1, return a -1
            if (leftHeight == -1)
                return -1;
            int rightHeight = height(root.right);
            if (rightHeight == -1)
                return -1;

            // if at any point, the difference between left and right subtree
            // height is more than 1, return a -1
            if (Math.abs(leftHeight - rightHeight) > 1)
                return -1;
            // or else return the height of the tree
            return Math.max(leftHeight, rightHeight) + 1;
        }
    }

    // ^ LCA in Binary Tree

    // ^ Check if two trees are identical or not

    // Given the roots of two binary trees p and q, write a function to check if
    // they are the same or not.

    // Two binary trees are considered the same if they are structurally identical,
    // and the nodes have the same value.

    class two_trees_are_identical_or_not {
        public boolean isSameTree(TreeNode p, TreeNode q) {

            if (p == null || q == null)
                return p == q;

            return (p.val == q.val) && isSameTree(p.left, q.left) && isSameTree(p.right, q.right);
        }
    }

    // ^ Zig Zag Traversal of Binary Tree

    // Given the root of a binary tree, return the zigzag level order traversal of
    // its nodes' values. (i.e., from left to right, then right to left for the next
    // level and alternate between).

    class zig_zag_traversal {
        public List<List<Integer>> zigzagLevelOrder(TreeNode root) {
            List<List<Integer>> ans = new LinkedList<>();
            if (root == null) {
                return ans;
            }
            Queue<TreeNode> queue = new LinkedList<>();
            queue.add(root);
            boolean odd = true;
            while (!queue.isEmpty()) {
                int size = queue.size();
                LinkedList<Integer> res = new LinkedList<>();
                for (int i = 1; i <= size; i++) {
                    // store poll in a variable, or else it will get null
                    // one poll value is used for the entire loop
                    // if we use poll again and again instead of node
                    // an element will get deleted each time
                    TreeNode node = queue.poll();

                    if (node.left != null) {
                        queue.add(node.left);
                    }
                    if (node.right != null) {
                        queue.add(node.right);
                    }
                    if (odd) {
                        res.add(node.val);
                    } else {
                        res.addFirst(node.val);
                    }

                }
                ans.add(res);
                odd = !odd;
            }
            return ans;
        }
    }

    // ^ Boundary Traversal of Binary Tree

    public class boundary_traversal {

        public boolean isLeaf(TreeNode node) {
            if (node.left == null && node.right == null)
                return true;
            else
                return false;
        }

        public ArrayList<Integer> traverseBoundary(TreeNode root) {
            ArrayList<Integer> list = new ArrayList<>();
            // if root node is not leaf, then we add it to the list
            if (isLeaf(root) == false)
                list.add(root.val);
            addLeft(root, list); // first ad all the left boundary elements
            addLeaves(root, list); // then add all the leaf nodes
            addRight(root, list); // then add all the right boundary elements

            return list;
        }

        public void addLeft(TreeNode root, ArrayList<Integer> res) {
            TreeNode curr = root.left;
            while (curr != null) {
                if (isLeaf(curr) == false)
                    res.add(curr.val);
                if (curr.left != null)
                    curr = curr.left;
                else
                    curr = curr.right;
            }
        }

        public void addLeaves(TreeNode root, ArrayList<Integer> res) {
            if (isLeaf(root) == true) {
                res.add(root.val);
                return;
            }
            if (root.left != null)
                addLeaves(root.left, res);
            if (root.right != null)
                addLeaves(root.right, res);
        }

        public void addRight(TreeNode root, ArrayList<Integer> res) {
            TreeNode curr = root.right;
            List<Integer> temp = new ArrayList<>();
            while (curr != null) {
                if (isLeaf(curr) == false)
                    temp.add(curr.val);
                if (curr.right != null)
                    curr = curr.right;
                else
                    curr = curr.left;
            }

            // for adding elements in reverse order
            for (int i = temp.size() - 1; i >= 0; --i) {
                res.add(temp.get(i));
            }
        }
    }

    // ^ Maximum path sum

    // ^ Construct Binary Tree from inorder and preorder

    // ^ Construct Binary Tree from Inorder and Postorder

    // ^ Symmetric Binary Tree

    // ^ Flatten Binary Tree to LinkedList

    // ^ Check if Binary Tree is the mirror of itself or not

    // ^ Check for Children Sum Property
}
