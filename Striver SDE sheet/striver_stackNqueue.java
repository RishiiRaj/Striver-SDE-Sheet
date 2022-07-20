//! Striver SDE sheet all Stack and Queue questions

import java.util.*;

class striver_stackNqueue {

    // ^ Implement Stack Using Arrays

    // ^ Implement Queue Using Arrays

    // ^ Implement Stack using Queue (using single queue)

    // Implement a last-in-first-out (LIFO) stack using only two queues. The
    // implemented stack should support all the functions of a normal stack (push,
    // top, pop, and empty).

    // Implement the MyStack class:

    // void push(int x) Pushes element x to the top of the stack.
    // int pop() Removes the element on the top of the stack and returns it.
    // int top() Returns the element on the top of the stack.
    // boolean empty() Returns true if the stack is empty, false otherwise.

    class MyStack {

        public MyStack() {

        }

        // use queue linked List and not normal queue
        private Queue<Integer> q = new LinkedList<>();

        public void push(int x) {
            // add an element normally
            q.add(x);
            // and remove all the rest of the elements and add them again
            for (int i = 1; i < q.size(); i++) {
                q.add(q.remove());
            }
        }

        public int pop() {
            return q.remove();
        }

        public int top() {
            return q.peek();
        }

        public boolean empty() {
            return q.isEmpty();
        }
    }

    /**
     * Your MyStack object will be instantiated and called as such:
     * MyStack obj = new MyStack();
     * obj.push(x);
     * int param_2 = obj.pop();
     * int param_3 = obj.top();
     * boolean param_4 = obj.empty();
     */

    // ^ Implement Queue using Stack (0(1) amortized method)

    // Implement a first in first out (FIFO) queue using only two stacks. The
    // implemented queue should support all the functions of a normal queue (push,
    // peek, pop, and empty).

    // Implement the MyQueue class:

    // void push(int x) Pushes element x to the back of the queue.
    // int pop() Removes the element from the front of the queue and returns it.
    // int peek() Returns the element at the front of the queue.
    // boolean empty() Returns true if the queue is empty, false otherwise.

    class MyQueue {

        // make 2 stacks, input and output

        private Stack<Integer> input = new Stack<>();
        private Stack<Integer> output = new Stack<>();

        public MyQueue() {

        }

        // push normally to input stack
        public void push(int x) {
            input.push(x);
        }

        // if output not empty then output.pop
        // else transfer all input stack elements to output and then putput.pop
        public int pop() {
            if (!output.isEmpty()) {
                return output.pop();
            } else {
                while (!input.isEmpty()) {
                    output.push(input.pop());
                }
                return output.pop();
            }
        }

        // if output not empty, return output top
        // else transfer all input elements to output, then output.top
        public int peek() {
            if (!output.isEmpty()) {
                return output.peek();
            } else {
                while (!input.isEmpty()) {
                    output.push(input.pop());
                }
                return output.peek();
            }
        }

        // normally check if both are empty
        public boolean empty() {
            return input.isEmpty() && output.isEmpty();
        }
    }

    /**
     * Your MyQueue object will be instantiated and called as such:
     * MyQueue obj = new MyQueue();
     * obj.push(x);
     * int param_2 = obj.pop();
     * int param_3 = obj.peek();
     * boolean param_4 = obj.empty();
     */

    // ^ Check for balanced parentheses

    // Given a string s containing just the characters '(', ')', '{', '}', '[' and
    // ']', determine if the input string is valid.

    // An input string is valid if:

    // Open brackets must be closed by the same type of brackets.
    // Open brackets must be closed in the correct order.

    // Example 1:

    // Input: s = "()"
    // Output: true

    class balanced_parenthesis {
        public boolean isValid(String str) {
            Stack<Character> s = new Stack<>();
            for (int i = 0; i < str.length(); i++) {
                char c = str.charAt(i);
                if (c == '(' || c == '{' || c == '[')
                    s.push(c);
                else if (s.isEmpty())
                    return false;
                else if (c == ')' && s.pop() != '(')
                    return false;
                else if (c == '}' && s.pop() != '{')
                    return false;
                else if (c == ']' && s.pop() != '[')
                    return false;
            }
            return s.isEmpty();
        }
    }

    // ^ Next Greater Element

    // The next greater element of some element x in an array is the first greater
    // element that is to the right of x in the same array.

    // You are given two distinct 0-indexed integer arrays nums1 and nums2, where
    // nums1 is a subset of nums2.

    // For each 0 <= i < nums1.length, find the index j such that nums1[i] ==
    // nums2[j] and determine the next greater element of nums2[j] in nums2. If
    // there is no next greater element, then the answer for this query is -1.

    // Return an array ans of length nums1.length such that ans[i] is the next
    // greater element as described above.

    // Example 1:

    // Input: nums1 = [4,1,2], nums2 = [1,3,4,2]
    // Output: [-1,3,-1]
    // Explanation: The next greater element for each value of nums1 is as follows:
    // - 4 is underlined in nums2 = [1,3,4,2]. There is no next greater element, so
    // the answer is -1.
    // - 1 is underlined in nums2 = [1,3,4,2]. The next greater element is 3.
    // - 2 is underlined in nums2 = [1,3,4,2]. There is no next greater element, so
    // the answer is -1.

    class next_greater_element {
        public int[] nextGreaterElement(int[] nums1, int[] nums2) {
            Map<Integer, Integer> map = new HashMap<>(); // map from x to next greater element of x
            Stack<Integer> stack = new Stack<>();
            for (int num : nums2) {
                while (!stack.isEmpty() && stack.peek() < num)
                    map.put(stack.pop(), num);
                stack.push(num);
            }
            for (int i = 0; i < nums1.length; i++)
                nums1[i] = map.getOrDefault(nums1[i], -1);
            return nums1;
        }
    }

    // ^ Sort a Stack

    // just use Collections.sort

    // ^ Next Smaller Element

    // ^ LRU cache (IMPORTANT)

    // ^ LFU Cache

    // ^ Largest rectangle in a histogram Link 1

    // ^ Sliding Window maximum

    // ^ Implement Min Stack

    // Design a stack that supports push, pop, top, and retrieving the minimum
    // element in constant time.

    // Implement the MinStack class:

    // MinStack() initializes the stack object.
    // void push(int val) pushes the element val onto the stack.
    // void pop() removes the element on the top of the stack.
    // int top() gets the top element of the stack.
    // int getMin() retrieves the minimum element in the stack.

    class min_stack {
        int min = Integer.MAX_VALUE;
        Stack<Integer> stack = new Stack<Integer>();

        public void push(int x) {
            // only push the old minimum value when the current
            // minimum value changes after pushing the new value x
            if (x <= min) {
                stack.push(min);
                min = x;
            }
            stack.push(x);
        }

        public void pop() {
            // if pop operation could result in the changing of the current minimum value,
            // pop twice and change the current minimum value to the last minimum value.
            if (stack.pop() == min)
                min = stack.pop();
        }

        public int top() {
            return stack.peek();
        }

        public int getMin() {
            return min;
        }
    }

    /**
     * Your MinStack object will be instantiated and called as such:
     * MinStack obj = new MinStack();
     * obj.push(val);
     * obj.pop();
     * int param_3 = obj.top();
     * int param_4 = obj.getMin();
     */

    // ^ Rotten Orange (Using BFS)

    // You are given an m x n grid where each cell can have one of three values:

    // 0 representing an empty cell,
    // 1 representing a fresh orange, or
    // 2 representing a rotten orange.
    // Every minute, any fresh orange that is 4-directionally adjacent to a rotten
    // orange becomes rotten.

    // Return the minimum number of minutes that must elapse until no cell has a
    // fresh orange. If this is impossible, return -1.

    // Example 1:

    // Input: grid = [[2,1,1],[1,1,0],[0,1,1]]
    // Output: 4

    class rotten_oranges {
        public int orangesRotting(int[][] grid) {
            if (grid == null || grid.length == 0)
                return 0;
            int rows = grid.length;
            int cols = grid[0].length;
            Queue<int[]> queue = new LinkedList<>();
            int oranges = 0;

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    // if rotten oranges found, add it to the queue
                    if (grid[i][j] == 2)
                        queue.offer(new int[] { i, j });
                    // just count the number of oranges
                    if (grid[i][j] != 0)
                        oranges++;
                }
            }

            if (oranges == 0)
                return 0; // if there are no oranges
            int minutes = 0, count = 0;
            int[] dx = { 0, 0, 1, -1 }; // x direction movement
            int[] dy = { 1, -1, 0, 0 }; // y direction movement

            // bfs starting from initially rotten oranges
            while (!queue.isEmpty()) {
                int size = queue.size(); // current size of queue
                count += size; // no. of rotten oranges according to queue

                for (int i = 0; i < size; i++) {
                    // for every queue element
                    int[] point = queue.poll();
                    for (int j = 0; j < 4; j++) // 4 because there are 4 direction to go in coordinate axis
                    {
                        int x = point[0] + dx[j];
                        int y = point[1] + dy[j];

                        // checking if x or y movement is not inside the grid or if the next cell is
                        // empty, or has an already rotten orange, then we skip the current iteration
                        if (x < 0 || y < 0 || x >= rows || y >= cols || grid[x][y] == 0 || grid[x][y] == 2)
                            continue;

                        grid[x][y] = 2; // make it as rotten
                        queue.offer(new int[] { x, y }); // and add the new rotten orange to the queue
                    }
                }
                // if we still have rotten oranges
                if (queue.size() != 0)
                    minutes++;
            }
            // if the original number of oranges and the oranges traversed by queue is same
            // it means all oranges have been rotten, else all oranges can't be rotten
            return oranges == count ? minutes : -1;
        }
    }

    // ^ Stock Span Problem

    // ^ Find the maximum of minimums of every window size

    // ^ The Celebrity Problem

}