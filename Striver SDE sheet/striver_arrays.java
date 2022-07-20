// ! Striver's SDE sheet all Arrays solutions

import java.util.*;

class striver_arrays {

    // ? Day 1 - Array's

    // ^ Set Matrix Zeroes
    // Given an m x n integer matrix matrix, if an element is 0, set its entire row
    // and column to 0's.
    // You must do it in place.

    class matrix_zeroes {
        public void setZeroes(int[][] matrix) {
            int rows = matrix.length, cols = matrix[0].length;

            // we take reference array's for row and column
            int d1[] = new int[rows];
            int d2[] = new int[cols];

            // fill the reference array's with -1
            Arrays.fill(d1, -1);
            Arrays.fill(d2, -1);

            // iterate through the matrix and where ever we find zero, we make both the
            // reference arrays as zero there
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (matrix[i][j] == 0) {
                        d1[i] = 0;
                        d2[j] = 0;
                    }
                }
            }
            // we check wherever the ref array is zero, we make the matrix element zero
            // there
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (d1[i] == 0 || d2[j] == 0) {
                        matrix[i][j] = 0;
                    }
                }
            }
        }
    }

    // ^ Pascal's triangle
    // Given an integer numRows, return the first numRows of Pascal's triangle.
    // In Pascal's triangle, each number is the sum of the two numbers directly
    // above it and number of rows is always equal to number of columns.

    class pascals_triangle {
        public List<List<Integer>> generate(int numRows) {
            // pascal's triangle is represented as a list of lists.
            List<List<Integer>> res = new ArrayList<List<Integer>>();
            List<Integer> row, pre = null;
            // row will be current row and pre will be previous row
            for (int i = 0; i < numRows; ++i) { // loop to iterate the number of time as the size of pascal's triangle
                // every time we visit a new row, we initialize it with an empty array list
                row = new ArrayList<Integer>();
                for (int j = 0; j <= i; ++j) { // number of columns will always be equal to number of rows, this is
                                               // always true in pascals triangle
                    // element is set to 1 if it is a border element
                    // since border elements are always 1 in a pascal's triangle
                    if (j == 0 || j == i)
                        row.add(1);
                    else // the sum of above 2 elements is added
                        row.add(pre.get(j - 1) + pre.get(j));
                }
                pre = row;
                // pre is made the current row, and new row will be made at the starting of the
                // loop
                res.add(row); // this list as added to the list of lists
            }
            return res;
        }
    }

    // ^ Next Permutation
    // A permutation of an array of integers is an arrangement of its members into a
    // sequence or linear order.
    // For example, for arr = [1,2,3], the following are considered permutations of
    // arr: [1,2,3], [1,3,2], [3,1,2], [2,3,1].

    // The next permutation of an array of integers is the next lexicographically
    // greater permutation of its integer. More formally, if all the permutations of
    // the array are sorted in one container according to their lexicographical
    // order, then the next permutation of that array is the permutation that
    // follows it in the sorted container. If such arrangement is not possible, the
    // array must be rearranged as the lowest possible order (i.e., sorted in
    // ascending order).

    class next_permutation {
        public void nextPermutation(int[] A) {

            // Algorithm-
            // linearly traverse from the back.
            // find first index i where A[i] < A[i+1]
            // again traverse from the back and find the first element which
            // is greater than the index i and store it in index j
            // swap values at index i and j
            // reverse from i+1 till last

            // base case - if array is empty or only one element is there
            if (A == null || A.length <= 1)
                return;
            int i = A.length - 2; // since we are traversing from the back, and we
            // will compare A[i] and A[i+1] so we start from length-2;
            // finding break point, i.e A[i]>A[i+1] traversing from the back
            while (i >= 0 && A[i] >= A[i + 1])
                i--; // i will be the break point.
            if (i >= 0) { // if there is a break point
                int j = A.length - 1;
                while (A[j] <= A[i])
                    j--; // we try to find someone who is greater then A[i]
                swap(A, i, j); // after finding it, we swap
            }
            // if there is no break point, the above block of code is never performed
            // and then reverse the right half, from the break point.
            reverse(A, i + 1, A.length - 1);
        }

        public void swap(int[] A, int i, int j) {
            int tmp = A[i];
            A[i] = A[j];
            A[j] = tmp;
        }

        public void reverse(int[] A, int i, int j) {
            while (i < j)
                swap(A, i++, j--);
        }
    }

    // ^ Maximum sum of sub-array ( Kadane's Algorithm )
    // Given an integer array nums, find the contiguous sub-array (containing at
    // least one number) which has the largest sum and return its sum.
    // A sub-array is a contiguous part of an array.

    class max_sum_subarray {
        public int maxSubArray(int[] nums) {
            int sum = 0;
            int max = Integer.MIN_VALUE;
            // we run a for each loop
            for (int i : nums) {
                sum += i; // store the current sum of the array
                if (sum > max)
                    max = sum; // store the max sum of the array
                if (sum < 0)
                    sum = 0; // if sum goes below zero, make sum as 0
            }
            return max;
        }
    }

    // ^ Sort colors OR Sort array of 0's 1's and 2's
    // Given an array nums with n objects colored red, white, or blue, sort them
    // in-place so that objects of the same color are adjacent, with the colors in
    // the order red, white, and blue.
    // We will use the integers 0, 1, and 2 to represent the color red, white, and
    // blue, respectively.
    // You must solve this problem without using the library's sort function.

    // basic approach, counting no. of zero's one's and tow's and replacing them

    class sort_colors_basic_approach {
        public void sortColors(int[] nums) {
            int zero = 0, one = 0, two = 0;
            for (int i = 0; i < nums.length; i++) {
                if (nums[i] == 0)
                    zero++;
                else if (nums[i] == 1)
                    one++;
                else
                    two++;
            }
            int i = 0;
            while (i < zero) {
                nums[i++] = 0;
            }
            while (i < one + zero) {
                nums[i++] = 1;
            }
            while (i < nums.length) {
                nums[i++] = 2;
            }
        }
    }

    // Dutch national flag algorithm

    class sort_colors {

        // IDK why but this swap function that I made, didn't work
        // because we have to pass the array as well inn the function parameter
        // public static void swap(int a, int b){
        // int temp=a;
        // a=b;
        // b=temp;
        // }
        public void sortColors(int[] nums) {
            int low = 0, mid = 0, high = nums.length - 1, temp;
            while (mid <= high) { // iterate mid
                switch (nums[mid]) {
                    case 0: { // if mid points to zero
                        // swap(nums[low], nums[mid]);
                        temp = nums[low];
                        nums[low] = nums[mid];
                        nums[mid] = temp;
                        low++; // move low one place ahead
                        mid++; // move mid also one place ahead
                        break;
                    }

                    case 1: { // if mid points to 1, don't do anything
                        mid++; // just move mid pointer ahead
                        break;
                    }

                    case 2: { // if mid points to 2
                        // swap(nums[mid], nums[high]);
                        temp = nums[mid];
                        nums[mid] = nums[high];
                        nums[high] = temp;
                        high--; // don't move the mid pointer
                        break;
                    }
                }
            }
        }
    }

    // ^ Stock buy and sell
    // You are given an array prices where prices[i] is the price of a given stock
    // on the ith day.
    // You want to maximize your profit by choosing a single day to buy one stock
    // and choosing a different day in the future to sell that stock.
    // Return the maximum profit you can achieve from this transaction. If you
    // cannot achieve any profit, return 0.

    public class stocks_buy_sell {
        // we traverse through the array, and keep a note of the smallest number
        // then as we traverse, we minus the min value from the current value
        // and similarly store the current max profit
        // all this is done in single traversal
        public int maxProfit(int prices[]) {
            int min_price = Integer.MAX_VALUE;
            int max_profit = 0;
            for (int i = 0; i < prices.length; i++) {
                if (prices[i] < min_price)
                    min_price = prices[i];
                else if (prices[i] - min_price > max_profit)
                    max_profit = prices[i] - min_price;
            }
            return max_profit;
        }
    }

    // ? Day 2 - Array's

    // ^ Rotate matrix 90 degree
    // You are given an n x n 2D matrix representing an image, rotate the image by
    // 90 degrees (clockwise).
    // You have to rotate the image in-place, which means you have to modify the
    // input 2D matrix directly. DO NOT allocate another 2D matrix and do the
    // rotation.
    // Input: matrix = [[1,2,3],[4,5,6],[7,8,9]]
    // Output: [[7,4,1],[8,5,2],[9,6,3]]

    class rotate_matrix {
        public void rotate(int[][] matrix) {
            // Algorithm is first to take transpose of the matrix and reverse every row

            // taking transpose
            for (int i = 0; i < matrix.length; i++) {
                for (int j = i; j < matrix[0].length; j++) {
                    int temp = matrix[i][j];
                    matrix[i][j] = matrix[j][i];
                    matrix[j][i] = temp;
                }
            }

            // take special care about the for loop conditions
            // since rows are now columns, we modify jth block of matrix, and not ith block
            // reversing every row using 2 pointer method
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length / 2; j++) { // j loop starts from 0, not from i and is till half of
                                                              // the matrix length
                    int temp = matrix[i][j];
                    matrix[i][j] = matrix[i][matrix.length - j - 1];
                    matrix[i][matrix.length - 1 - j] = temp;
                }
            }
        }
    }

    // ^ Merge Overlapping Intervals
    // Given an array of intervals where intervals[i] = [starti, endi], merge all
    // overlapping intervals, and return an array of the non-overlapping intervals
    // that cover all the intervals in the input.
    // Input: intervals = [[1,3],[2,6],[8,10],[15,18]]
    // Output: [[1,6],[8,10],[15,18]]
    // Explanation: Since intervals [1,3] and [2,6] overlaps, merge them into [1,6].

    class merge_overlapping_intervals {
        public int[][] merge(int[][] intervals) {

            // base case - if given array is null
            if (intervals == null || intervals.length == 0)
                return intervals;

            // sort intervals by starting value
            Arrays.sort(intervals, (a, b) -> Integer.compare(a[0], b[0]));

            // if end of previous interval is more than the start of current interval then
            // there is a overlap
            LinkedList<int[]> mergedIntervals = new LinkedList<>();
            for (int[] curr : intervals) {
                // if list empty or no overlap simply add current interval
                // if end of an interval is less than start of next element, then there is no
                // overlap
                if (mergedIntervals.isEmpty() || mergedIntervals.getLast()[1] < curr[0])
                    mergedIntervals.add(curr);
                // else if overlap exists then merge current interval with the previous interval
                else
                    // we modify the last element in the list and make it an overlapping interval
                    mergedIntervals.getLast()[1] = Math.max(mergedIntervals.getLast()[1], curr[1]);
            }
            return mergedIntervals.toArray(new int[0][]);
        }
    }

    // ^ Merge 2 Sorted Arrays without Extra Space
    // You are given two integer arrays nums1 and nums2, sorted in non-decreasing
    // order, and two integers m and n, representing the number of elements in nums1
    // and nums2 respectively.
    // Merge nums1 and nums2 into a single array sorted in non-decreasing order.
    // The final sorted array should not be returned by the function, but instead be
    // stored inside the array nums1. To accommodate this, nums1 has a length of m +
    // n, where the first m elements denote the elements that should be merged, and
    // the last n elements are set to 0 and should be ignored. nums2 has a length of
    // n.

    class merge_without_extra_space {
        public void merge(int[] nums1, int m, int[] nums2, int n) {
            // here m is not the size of first array, but the number of elements
            // in first array, size of 2nd array is m+n
            int i, j, k;
            if (m != 0 && n != 0) { // if none of the arrays are empty
                for (i = 0; i < m; i++) {
                    // if ith element in first array is greater than first element in 2nd array
                    // we swap them and then sort the second array.
                    if (nums1[i] > nums2[0]) {
                        int temp = nums1[i];
                        nums1[i] = nums2[0];
                        nums2[0] = temp;
                    }
                    // insertion sort is used here
                    // we can use Arrays.sort here as well.
                    int first = nums2[0];
                    for (k = 1; k < n && nums2[k] < first; k++) {
                        nums2[k - 1] = nums2[k];
                    }
                    nums2[k - 1] = first;
                }
                // after sorting, we copy the 2nd array in the first array in the
                // remaining spaces in array1 which were filled with zero.
                for (j = m, i = 0; j < m + n && i < n; j++, i++) {
                    nums1[j] = nums2[i];
                }
            }
            // if any one array is empty
            else {
                for (j = 0, i = 0; j < n; j++, i++) {
                    nums1[i] = nums2[j];
                }
            }
        }
    }

    // ^ Find Duplicate number in array of N+1 integers
    // Given an array of integers nums containing n + 1 integers where each integer
    // is in the range [1, n] inclusive.
    // There is only one repeated number in nums, return this repeated number.
    // You must solve the problem without modifying the array nums and uses only
    // constant extra space.

    // Naive approach is to sort and check for consecutive elements to be same or
    // not

    // better approach is to use HashSet.

    class find_duplicate {
        // we use the tortoise and hare method.
        public int findDuplicate(int[] nums) {
            int slow = nums[0], fast = nums[0];
            do {
                slow = nums[slow]; // move slow pointer one step
                fast = nums[nums[fast]]; // move fast pointer 2 steps
            } while (slow != fast); // until they meet.

            fast = nums[0]; // initialize fast to beginning of array.

            while (fast != slow) {
                // increment both by one step until they meet
                slow = nums[slow];
                fast = nums[fast];
            }
            // return any one, as that is the answer.
            return slow;
        }
    }

    // ^ Repeat and missing number
    // You are given an array of size ‘N’. The elements of the array are in the
    // range from 1 to ‘N’.
    // Ideally, the array should contain elements from 1 to ‘N’. But due to some
    // miscalculations, there is a number R in the range [1, N] which appears in the
    // array twice and another number M in the range [1, N] which is missing from
    // the array.
    // Your task is to find the missing number (M) and the repeating number (R).
    // Consider an array of size six. The elements of the array are
    // { 6, 4, 3, 5, 5, 1 }.
    // The array should contain elements from one to six. Here, 2 is not present and
    // 5 is occurring twice. Thus, 2 is the missing number (M) and 5 is the
    // repeating number (R).

    // Naive approach is to set up a frequency array
    // after updating the frequency array, we traverse the freq. array, and the
    // index which has value 2 is the repeating number, and the index having value 0
    // other then 0th index is the missing number.

    // optimized approach is to use XOR property.

    // ^ Inversion of an array
    // For a given integer array/list 'ARR' of size 'N', find the total number of
    // 'Inversions' that may exist.
    // An inversion is defined for a pair of integers in the array/list when the
    // following two conditions are met.
    // A pair ('ARR[i]', 'ARR[j]') is said to be an inversion when:
    // 1. 'ARR[i] > 'ARR[j]'
    // 2. 'i' < 'j'
    // Where 'i' and 'j' denote the indices ranging from [0, 'N').
    // Sample Input 1 :
    // 3
    // 3 2 1
    // Sample Output 1 :
    // 3
    // Explanation Of Sample Output 1:
    // We have a total of 3 pairs which satisfy the condition of inversion. (3, 2),
    // (2, 1) and (3, 1).

    // ? Day 3 - Arrays

    // ^ Search in a 2D matrix (LeetCode)
    // Write an efficient algorithm that searches for a value target in an m x n
    // integer matrix matrix. This matrix has the following properties:
    // Integers in each row are sorted from left to right.
    // The first integer of each row is greater than the last integer of the
    // previous row.

    // Naive approach is to do linear search. O(N^2)
    // better approach is to do binary search in each row since rows are sorted.
    // O(NlogM)

    class search_in_matrix {
        public boolean searchMatrix(int[][] matrix, int target) {
            // since rows are sorted
            int i = 0, j = matrix[0].length - 1, n = matrix.length;
            // we start by placing the pointer to the last index of the first row
            // everything at the left of this index will be lesser than that value
            // everything below that index in that column, will be greater than that index
            // value
            while (i < n && j >= 0) {
                if (matrix[i][j] == target) {
                    return true;
                }
                // if target is less than that the index value, we move one column left
                if (matrix[i][j] > target) {
                    j--;
                }
                // if target is greater than the index value, then we are sure that it is in
                // that column only, then we need to increment the row
                else {
                    i++;
                }
            }
            // if target is not in the matrix, pointer will keep changing and eventually
            // get out of the while loop.
            return false;
        }
    }

    // most optimal approach is to use binary search, since the array is sorted.
    class search_in_matrix_optimal {
        public boolean searchMatrix(int[][] matrix, int target) {
            if (matrix.length == 0)
                return false;

            // if we take values of the in an order, it will be a sorted array
            // we will perform binary search in that array without using extra space
            // just by using the index values.

            int n = matrix.length;
            int m = matrix[0].length;

            int low = 0, high = (m * n) - 1;

            // just simple binary search;
            while (low <= high) {
                int mid = low + (high - low) / 2;
                // finding value of row and column by
                // row value = mid value / column length
                // column value = mid value % column length
                int mid_row = mid / m;
                int mid_col = mid % m;

                if (matrix[mid_row][mid_col] == target)
                    return true;

                if (matrix[mid_row][mid_col] < target)
                    low = mid + 1;

                else
                    high = mid - 1;
            }
            return false;
        }
    }

    // ^ Pow(x,n)
    // Implement pow(x, n), which calculates x raised to the power n (i.e., xn).
    // Input: x = 2.00000, n = 10
    // Output: 1024.00000
    // Input: x = 2.00000, n = -2
    // Output: 0.25000
    // Explanation: 2-2 = 1/22 = 1/4 = 0.25

    // brute force approach uses O(n) time complexity, gives TLE
    // this approach is logN
    class power_x_to_n {
        public double myPow(double x, int n) {
            double ans = 1.0;
            // storing n in a long type to escape an edge case
            // edge case is that int ranges from -2,147,483,648 to 2,147,483,647
            // last number cannot be converted to int, it will overflow
            // this is the edge case, and we need to take long datatype
            long nn = n;
            // if n is negative, we convert it to positive, and do the calculation;
            if (nn < 0)
                nn = -1 * nn;

            // we will keep decreasing the value of n until it is zero
            while (nn > 0) {
                // if nn is an odd power then multiply x with ans ans reduce nn by 1
                if (nn % 2 == 1) {
                    ans = ans * x;
                    nn = nn - 1;
                } else { // Else multiply x with itself and divide nn by two.
                    x = x * x;
                    nn = nn / 2;
                }
            }
            if (n < 0)
                ans = (double) (1.0) / (double) (ans);
            return ans;
        }
    }

    // ^ Majority Element (>N/2 Times)
    // Given an array nums of size n, return the majority element.
    // The majority element is the element that appears more than ⌊n / 2⌋ times. You
    // may assume that the majority element always exists in the array.

    // brute force solution
    class majority_elem_half {
        public int count(int arr[], int x) {
            int c = 0;
            for (int i = 0; i < arr.length; i++) {
                if (arr[i] == x)
                    c++;
            }
            return c;
        }

        public int majorityElement(int[] nums) {

            int l = nums.length;
            for (int i = 0; i < l; i++) {
                if (count(nums, nums[i]) > l / 2) {
                    return nums[i];
                }

            }
            return l;
        }
    }

    // optimized solution could be using HashMap or frequency array

    // Optimal Approach
    class majority_elem_half_optimized {
        public int majorityElement(int[] nums) {
            // we use moore's voting algorithm
            // majority element gets cancelled by the minority element
            // we use 2 variables, one for storing an element, and one for its frequency
            // the last element remaining in the variable, will be the majority element

            int count = 0, majority_element = 0;
            // linearly traversing in the array
            for (int num : nums) {
                // if count reaches zero, we set the current array element as majority element
                if (count == 0)
                    majority_element = num;
                if (num == majority_element)
                    count++; // count its occurrence
                else
                    count--; // if the current element does not match, count is decreased by one
                // this is how the majority elements gets cancelled out by the minority
                // elements.
            }
            return majority_element;
        }
    }

    // ^ Majority Element (>N/3 Times)
    // Given an integer array of size n, find all elements that
    // appear more than n/3 times
    // Example 1:
    // Input: nums = [3,2,3]
    // Output: [3]

    // brute force approach is to count occurrence to each element and which all
    // have appeared more than 1/3 times, we add them to the resultant list.

    // optimized approach is to use HashMap
    class majority_elem_one_third_optimized {
        public List<Integer> majorityElement(int[] arr) {
            HashMap<Integer, Integer> mp = new HashMap<>();
            ArrayList<Integer> ans = new ArrayList<>();
            int n = arr.length;
            // add each element and its frequency in hash map
            for (int i = 0; i < n; i++) {
                mp.put(arr[i], mp.getOrDefault(arr[i], 0) + 1);
            }
            // checking which element is present more then 1/3 times
            // andadding it to the resultant list
            for (int x : mp.keySet()) {
                if (mp.get(x) > (n / 3))
                    ans.add(x);
            }

            return ans;
        }
    }

    // optimal approach is to use Moore's voting algorithm.
    class majority_elem_one_third {
        public List<Integer> majorityElement(int[] nums) {
            // moore's voting algorithm
            // at most, 2 majority element could be present, and at least zero
            int number1 = -1, number2 = -1, count1 = 0, count2 = 0, len = nums.length;
            for (int num : nums) {
                if (num == number1)
                    count1++;
                else if (num == number2)
                    count2++;
                else if (count1 == 0) {
                    number1 = num;
                    count1 = 1;
                } else if (count2 == 0) {
                    number2 = num;
                    count2 = 1;
                } else {
                    count1--;
                    count2--;
                }
            }
            // till here, elem1 and elem2 store the majority elements
            ArrayList<Integer> ans = new ArrayList<Integer>();
            count1 = 0;
            count2 = 0;
            for (int numz : nums) {
                if (numz == number1)
                    count1++;
                else if (numz == number2)
                    count2++;
            }
            // we again counted the frequency of elem1 and elem 2
            // check if its freq is more than l/3
            // if yes we add it to the list
            if (count1 > len / 3)
                ans.add(number1);
            if (count2 > len / 3)
                ans.add(number2);
            return ans;
        }
    }

    // ^ Grid unique paths
    // There is a robot on an m x n grid. The robot is initially located at the
    // top-left corner (i.e., grid[0][0]). The robot tries to move to the
    // bottom-right corner (i.e., grid[m - 1][n - 1]). The robot can only move
    // either down or right at any point in time.

    // Given the two integers m and n, return the number of possible unique paths
    // that the robot can take to reach the bottom-right corner.

    // The test cases are generated so that the answer will be less than or equal to
    // 2 * 109.

    // Example 2:

    // Input: m = 3, n = 2
    // Output: 3
    // Explanation: From the top-left corner, there are a total of 3 ways to reach
    // the bottom-right corner:
    // 1. Right -> Down -> Down
    // 2. Down -> Down -> Right
    // 3. Down -> Right -> Down

    // basic approach is to use recursion, and it uses exponential time and space
    class grid_unique_paths_recursion {
        public int uniquePaths(int m, int n) {
            // base case
            if (m == 1 || n == 1)
                return 1;

            // move down
            int downMove = uniquePaths(m - 1, n);
            // move right
            int rightMove = uniquePaths(m, n - 1);

            return downMove + rightMove;
        }
    }

    // better approach is to use dynamic programming and HashMap.
    class grid_unique_paths_better {
        private Map<String, Integer> map = new HashMap<String, Integer>();

        public int uniquePaths(int m, int n) {
            // base case
            if (m == 1 || n == 1)
                return 1;

            // check if we have already calculated unique paths for cell(m, n)
            String cell = new String(m + "," + n);
            // if yes, then get its value from our memoization table
            if (map.containsKey(cell))
                return map.get(cell);

            // else, explore the down move
            int downMove = uniquePaths(m - 1, n);
            // explore the right move
            int rightMove = uniquePaths(m, n - 1);

            // put the value obtained for unique paths from cell(m, n)
            map.put(cell, downMove + rightMove);

            return downMove + rightMove;
        }
    }
    // most optimal approach is to use combination

    // on observing, we notice that no. or right moves is one less than the number
    // of rows
    // and number of down moves is one less than the number of columns
    // total directions to choose is (m-1)+(n-1) = m+n-2
    // answer will be (m+n-2)C(n-1) or (m+n-2)C(m-1) both will give the same answer
    // shortcut to find out nCr is, find r! and number of multiplications in
    // denominator
    // that will be equal to number of multiplications is n! from the starting
    // example 10C3 = (10x9x8)/(3x2x1)

    class grid_unique_paths_optimal {
        public int uniquePaths(int m, int n) {
            if (m == 1 || n == 1)
                return 1;
            // answer will be (m+n-2)C(n-1) or (m+n-2)C(m-1)
            int nc = m + n - 2;
            int r = m - 1;
            double res = 1; // taking integer here will give wrong answer, idk why
            for (int i = 1; i <= r; i++) {
                // here i will be denominator
                // and nc-r+i will be numerator
                res = res * (nc - r + i) / i;
            }
            return (int) res;
        }
    }

    // ^ Reverse Pairs
    // Given an integer array nums, return the number of reverse pairs in the array.

    // A reverse pair is a pair (i, j) where 0 <= i < j < nums.length and nums[i] >
    // 2 * nums[j].

    // Example 1:

    // Input: nums = [1,3,2,3,1]
    // Output: 2

    // brute force approach is to use nested loops, and check for every pair
    // if the pair satisfies the condition, we increase counter by one.
    // but... hah! this is leet code hard, and will give TLE
    class reverse_pairs_brute {
        public int reversePairs(int[] arr) {
            int count = 0;
            for (int i = 0; i < arr.length; i++)
                for (int j = i + 1; j < arr.length; j++)
                    if ((long) arr[i] > (long) 2 * arr[j]) // 21,47,483,647, which is max value of int, was given in
                                                           // some test cases
                        count++;
            return count;
        }
    }

    // better approach is to use modified MergeSort.
    // didn't understand that so will do afterwards.

    // ^ 2 Sum problem
    // Given an array of integers nums and an integer target, return indices of the
    // two numbers such that they add up to target.

    // You may assume that each input would have exactly one solution, and you may
    // not use the same element twice.

    // You can return the answer in any order.

    // Example 1:

    // Input: nums = [2,7,11,15], target = 9
    // Output: [0,1]
    // Explanation: Because nums[0] + nums[1] == 9, we return [0, 1].

    // Brute force approach is to use nested loops and check for all possible pairs
    class two_sum_brute {
        public int[] twoSum(int[] nums, int target) {
            int ans[] = new int[2];
            for (int i = 0; i < nums.length; i++) {
                for (int j = i + 1; j < nums.length; j++) {
                    if (nums[i] + nums[j] == target) {
                        ans[0] = i;
                        ans[1] = j;
                    }
                }
            }
            return ans;
        }
    }

    // Optimized approach is to use HashMap
    class two_sum_optimal {
        public int[] twoSum(int[] nums, int target) {
            int[] ans = new int[2]; // make an array of size 2 to store the answer
            // create a HashMap for calculation
            Map<Integer, Integer> map = new HashMap<Integer, Integer>();
            int l = nums.length;
            // iterate through the array
            for (int i = 0; i < l; i++) {
                // if target-current array element is present in the hashmap, then those 2 are
                // the ans
                if (map.containsKey(target - nums[i])) {
                    ans[0] = i;
                    ans[1] = map.get(target - nums[i]);
                }
                // if not, we put the current array value in the hashmap
                map.put(nums[i], i);
            }
            return ans;
        }
    }

    // ^ 4 Sum problem
    // Given an array nums of n integers, return an array of all the unique
    // quadruplets [nums[a], nums[b], nums[c], nums[d]] such that:

    // 0 <= a, b, c, d < n
    // a, b, c, and d are distinct.
    // nums[a] + nums[b] + nums[c] + nums[d] == target
    // You may return the answer in any order.

    // Example 1:

    // Input: nums = [1,0,-1,0,-2,2], target = 0
    // Output: [[-2,-1,1,2],[-2,0,0,2],[-1,0,0,1]]

    // brute force approach is to sort the array
    // use 3 nested loops and use binary search for searching the 4th element
    // all possible solutions are put in a HashSet to find all the unique solutions

    // better approach is to first fort the arrya
    // then use 2 nested loops and use 2 pointers to find out the rest of the 2
    // elements
    class four_sum_optimal {
        public List<List<Integer>> fourSum(int[] nums, int target) {
            ArrayList<List<Integer>> res = new ArrayList<List<Integer>>();
            if (nums == null || nums.length == 0)
                return res;
            int l = nums.length;
            Arrays.sort(nums);
            for (int i = 0; i < l; i++) {
                // target minus the value at i pointer
                int target_3 = target - nums[i];
                for (int j = i + 1; j < l; j++) {
                    // target minus the value at i and j pointer
                    int target_2 = target_3 - nums[j];
                    // creating the 2 pointers to find the target_2
                    int front = j + 1;
                    int back = l - 1;
                    // iterate in the remaining elements
                    while (front < back) {
                        int two_sum = nums[front] + nums[back];
                        // if the target sum is greater than the sum we got, then in order to increase
                        // the sum, we increment the front pointer to get a greater value
                        // because the array is sorted
                        if (two_sum < target_2)
                            front++;
                        // if target sum is smaller, we move the back pointer to the left
                        else if (two_sum > target_2)
                            back--;
                        // else when the sum is satisfied (target and two_sum is equal)
                        else {
                            List<Integer> quad = new ArrayList<>();
                            quad.add(nums[i]);// quad posi 0
                            quad.add(nums[j]);// quad posi 1
                            quad.add(nums[front]);// quad posi 2
                            quad.add(nums[back]);// quad posi 3
                            res.add(quad);

                            // skipping duplicates for number 3
                            while (front < back && nums[front] == quad.get(2))
                                ++front;

                            // skipping duplicates fro number 4
                            while (front < back && nums[back] == quad.get(3))
                                --back;
                        }
                    }
                    // skipping duplicates for number 2
                    // j+1 taaki end tak na pahuch jaaye
                    // nums[j+1]==nums[j] issliye kyuki ek posi peeche tak rakhna hai, kyuki j wale
                    // for loop me ek increment hoga hi
                    while (j + 1 < l && nums[j + 1] == nums[j])
                        ++j;
                }
                // skipping duplicates for number 1
                while (i + 1 < l && nums[i + 1] == nums[i])
                    ++i;
            }
            return res;
        }
    }

    // ^ Longest Consecutive subsequence
    // Given an unsorted array of integers nums, return the length of the longest
    // consecutive elements sequence.

    // You must write an algorithm that runs in O(n) time.

    // Example 1:

    // Input: nums = [100,4,200,1,3,2]
    // Output: 4
    // Explanation: The longest consecutive elements sequence is [1, 2, 3, 4].
    // Therefore its length is 4.

    class Longest_consecutive_subsequence_brute {
        public int longestConsecutive(int[] nums) {
            if (nums.length == 0 || nums == null) {
                return 0;
            }

            // first we sort the array
            Arrays.sort(nums);

            int ans = 1; // initialize the current answer as one
            int prev = nums[0];
            int cur = 1;

            for (int i = 1; i < nums.length; i++) {
                if (nums[i] == prev + 1) { // checking if ith element is one more than the previous element
                    cur++; // we increment the current sequence length.
                } else if (nums[i] != prev) { // as the streak is broken, i.e. there is no consecutive number, current
                                              // is set
                                              // to one again
                    cur = 1;
                }
                prev = nums[i]; // prev is updated
                ans = Math.max(ans, cur); // store the maximum streak in ans comparing it with current streak
            }
            return ans;
        }
    }

    // the above discussed solution is the naive approach, the brute force solution

    // below solution is the optimized solution

    class Longest_consecutive_subsequence_optimal {
        public int longestConsecutive(int[] nums) {
            HashSet<Integer> h = new HashSet<>();
            // store all the elements in a hash set
            for (int i : nums) {
                h.add(i);
            }

            int longestStreak = 0;

            for (int i : nums) {
                if (!h.contains(i - 1)) { // if HashSet contains number one less than the current number then we iterate
                                          // forward, if it doesn't contain, then we enter this part of the code
                    int curr = i;
                    int currStreak = 1; // current streak length is always one

                    // we go till the smallest consecutive number and check for one one greater
                    // number is present or not, like this we will have to check for least number of
                    // times.
                    while (h.contains(curr + 1)) { // if HashSet contains one more than the current number
                        curr += 1; // we will next pass the next number, so one two+current number is calculated
                        currStreak += 1; // streak is incremented by one
                    }
                    longestStreak = Math.max(longestStreak, currStreak); // longest streak is maintains
                }
            }
            return longestStreak;
        }
    }

    // ^ Largest Sub-array with 0 Sum
    // Given an array having both positive and negative integers. The task is to
    // compute the length of the largest sub-array with sum 0.
    // Example 1:

    // Input:
    // N = 8
    // A[] = {15,-2,2,-8,1,7,10,23}
    // Output: 5
    // Explanation: The largest sub-array with
    // sum 0 will be -2 2 -8 1 7.

    // brute force approach is to find all sub-arrays
    // optimal approach is to use HashMap

    // ^ Count number of Sub-arrays with given XOR K

    // ^ Longest Substring without repeat

    // ^ 3 Sum

    // ^ Trapping Rainwater

    // ^ Remove duplicates from a sorted Array

    // Given an integer array nums sorted in non-decreasing order, remove the
    // duplicates in-place such that each unique element appears only once. The
    // relative order of the elements should be kept the same.

    // Since it is impossible to change the length of the array in some languages,
    // you must instead have the result be placed in the first part of the array
    // nums. More formally, if there are k elements after removing the duplicates,
    // then the first k elements of nums should hold the final result. It does not
    // matter what you leave beyond the first k elements.

    // Return k after placing the final result in the first k slots of nums.

    // Do not allocate extra space for another array. You must do this by modifying
    // the input array in-place with O(1) extra memory.

    // Example 1:

    // Input: nums = [1,1,2]
    // Output: 2, nums = [1,2,_]
    // Explanation: Your function should return k = 2, with the first two elements
    // of nums being 1 and 2 respectively.
    // It does not matter what you leave beyond the returned k (hence they are
    // underscores).

    class remove_duplicates_from_sorted_array {
        public int removeDuplicates(int[] nums) {
            if (nums.length == 0 || nums.length == 1) {
                return 1;
            }
            int i = 0; // i pinter at zero'th index
            // j pointer starts from 1st index
            for (int j = 1; j < nums.length; j++) {
                // if values are not equal, i pointer increments
                // and i index value is updated to j index value
                if (nums[i] != nums[j]) {
                    i++;
                    nums[i] = nums[j];
                }
            }
            // size of unique elements array is displayed
            return i + 1;
        }
    }

    // ^ Max Consecutive Ones

    // Given a binary array nums, return the maximum number of consecutive 1's in
    // the array.

    // Example 1:

    // Input: nums = [1,1,0,1,1,1]
    // Output: 3
    // Explanation: The first two digits or the last three digits are consecutive
    // 1s. The maximum number of consecutive 1s is 3.

    class max_consecutive_ones {
        public int findMaxConsecutiveOnes(int[] nums) {
            int c = 0, max = 0;
            for (int i = 0; i < nums.length; i++) {
                if (nums[i] == 1) {
                    c++;
                    max = Math.max(max, c);
                } else
                    c = 0;
            }
            return max;
        }
    }
}