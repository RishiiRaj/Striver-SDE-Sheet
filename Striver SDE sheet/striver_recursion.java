//! Striver SDE sheet, all recursion and backtracking solutions.

import java.util.*;

public class striver_recursion {

    // ? Day 9 recursion and backtracking

    // ^ Subset sums

    // Given a list arr of N integers, print sums of all subsets in it.

    // Example 1:
    // Input:
    // N = 2
    // arr[] = {2, 3}
    // Output:
    // 0 2 3 5
    // Explanation:
    // When no elements is taken then Sum = 0.
    // When only 2 is taken then Sum = 2.
    // When only 3 is taken then Sum = 3.
    // When element 2 and 3 are taken then
    // Sum = 2+3 = 5.

    class subset_sum {
        ArrayList<Integer> subsetSums(ArrayList<Integer> arr, int N) {
            ArrayList<Integer> res = new ArrayList<>();

            func(0, 0, arr, N, res);
            Collections.sort(res);
            return res;
        }

        public void func(int index, int sum, ArrayList<Integer> arr, int N, ArrayList<Integer> res) {
            if (index == N) {
                res.add(sum);
                return;
            }

            // pick up an element
            // index is incremented and sum is added
            func(index + 1, sum + arr.get(index), arr, N, res);
            // do not pick up an element
            // index is incremented but sum is not added
            func(index + 1, sum, arr, N, res);
        }
    }

    // ^ Subset 2

    // ^ Combination sum-1

    // Given an array of distinct integers candidates and a target integer target,
    // return a list of all unique combinations of candidates where the chosen
    // numbers sum to target. You may return the combinations in any order.

    // The same number may be chosen from candidates an unlimited number of times.
    // Two combinations are unique if the frequency of at least one of the chosen
    // numbers is different.

    // It is guaranteed that the number of unique combinations that sum up to target
    // is less than 150 combinations for the given input.

    // Example 1:

    // Input: candidates = [2,3,6,7], target = 7
    // Output: [[2,2,3],[7]]
    // Explanation:
    // 2 and 3 are candidates, and 2 + 2 + 3 = 7. Note that 2 can be used multiple
    // times.
    // 7 is a candidate, and 7 = 7.
    // These are the only two combinations.

    class combination_sum_1 {
        public List<List<Integer>> combinationSum(int[] candidates, int target) {
            List<List<Integer>> ans = new ArrayList<>();
            findCombi(0, candidates, target, ans, new ArrayList<Integer>());
            return ans;
        }

        public void findCombi(int index, int[] candidates, int target, List<List<Integer>> ans, List<Integer> ds) {
            // if we've reached the end, it means we have checked every index
            if (index == candidates.length) {
                // we keep reducing target by subtracting
                // once it reaches zero, it means we have found the correct combination
                // and we add that list to out answer list
                if (target == 0) {
                    ans.add(new ArrayList<>(ds));
                }
                return;
            }
            // if we decide to pick an element
            if (candidates[index] <= target) {
                ds.add(candidates[index]);
                // we do not move index, because that element can be re used
                findCombi(index, candidates, target - candidates[index], ans, ds);
                // we remove the temp addition in the DS, because it is a reference DS
                ds.remove(ds.size() - 1);
            }
            // if we decide to not pick an element
            findCombi(index + 1, candidates, target, ans, ds);
        }
    }

    // ^ Combination sum-2

    // Given a collection of candidate numbers (candidates) and a target number
    // (target), find all unique combinations in candidates where the candidate
    // numbers sum to target.

    // Each number in candidates may only be used once in the combination.

    // Note: The solution set must not contain duplicate combinations.

    // Example 1:

    // Input: candidates = [10,1,2,7,6,1,5], target = 8
    // Output:
    // [
    // [1,1,6],
    // [1,2,5],
    // [1,7],
    // [2,6]
    // ]

    class combination_sum_2 {
        public List<List<Integer>> combinationSum2(int[] candidates, int target) {
            Arrays.sort(candidates);
            List<List<Integer>> ans = new ArrayList<>();
            findCombinations2(0, candidates, target, ans, new ArrayList<>());
            return ans;
        }

        private void findCombinations2(int index, int[] arr, int target, List<List<Integer>> ans, List<Integer> temp) {
            if (target == 0) {
                ans.add(new ArrayList<>(temp));
                return;
            }

            for (int i = index; i < arr.length; i++) {
                // if i is the first element, it will be picked, irrespective of weather it is
                // a duplicate element ot not
                // thats why i>index check is done and next check is for duplicates
                if (i > index && arr[i] == arr[i - 1])
                    continue;
                // since array is sorted, the next elements will also be greater than rhe
                // target, so we break out
                if (arr[i] > target)
                    break;

                temp.add(arr[i]);
                findCombinations2(i + 1, arr, target - arr[i], ans, temp);
                temp.remove(temp.size() - 1);
            }
        }
    }

    // ^ Palindrome Partitioning

    // Given a string s, partition s such that every substring of the partition is a
    // palindrome. Return all possible palindrome partitioning of s.
    // A palindrome string is a string that reads the same backward as forward.

    // Example 1:

    // Input: s = "aab"
    // Output: [["a","a","b"],["aa","b"]]

    class palindrome_partitioning {
        public List<List<String>> partition(String s) {
            List<List<String>> ans = new ArrayList<>();
            List<String> temp = new ArrayList<>();
            backtrack(s, 0, temp, ans);
            return ans;
        }

        public void backtrack(String s, int index, List<String> temp, List<List<String>> ans) {
            // if the string has been partitioned completely, that's the base case
            if (index == s.length()) {
                ans.add(new ArrayList<>(temp));
                return;
            }

            // loop from index to end of string
            for (int i = index; i < s.length(); ++i) {
                // if it is a palindrome, we add it to list, and backtrack
                if (isPal(s, index, i)) {
                    temp.add(s.substring(index, i + 1));
                    backtrack(s, i + 1, temp, ans);
                    // and remove the last element at the end
                    temp.remove(temp.size() - 1);
                }
            }
        }

        public boolean isPal(String s, int start, int end) {
            while (start <= end) {
                if (s.charAt(start++) != s.charAt(end--))
                    return false;
            }
            return true;
        }
    }

    // ^ K-th permutation Sequence (LeetCode hard)

    // The set [1, 2, 3, ..., n] contains a total of n! unique permutations.

    // By listing and labeling all of the permutations in order, we get the
    // following sequence for n = 3:

    // "123"
    // "132"
    // "213"
    // "231"
    // "312"
    // "321"
    // Given n and k, return the kth permutation sequence.

    // Example 1:

    // Input: n = 3, k = 3
    // Output: "213"

    // this is optimal solution using mathematics, brute force is using recursion,
    // and finding all the possible permutations

    class kth_permutation {
        public String getPermutation(int n, int k) {
            int fact = 1;
            List<Integer> list = new ArrayList<>();
            // compute n-1 factorial, since loop is < n and not <= n
            for (int i = 1; i < n; i++) {
                fact = fact * i;
                list.add(i); // adding n-1 numbers in the list
            }
            list.add(n); // adding n'th number in the list.
            String ans = "";
            k = k - 1; // since it is zero based indexing, if we're asked to find the 17th permutation,
                       // we will find the 16th permutation
            // infinite loop
            while (true) {
                ans = ans + list.get(k / fact);
                list.remove(k / fact);
                if (list.size() == 0)
                    break; // till the list becomes zero;

                k = k % fact; // next value of k will be modulo
                fact = fact / list.size(); // get one factorial less, reducing factorial aat every step
            }
            return ans;
        }
    }

    // ^ Print all permutations of a string/array

    // Given an array nums of distinct integers, return all the possible
    // permutations. You can return the answer in any order.

    // Example 1:

    // Input: nums = [1,2,3]
    // Output: [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]]

    class permutations {
        public List<List<Integer>> permute(int[] nums) {
            List<List<Integer>> res = new ArrayList<>();
            List<Integer> ds = new ArrayList<>();
            boolean[] freq = new boolean[nums.length];
            backtrack(nums, res, ds, freq);
            return res;
        }

        private void backtrack(int[] nums, List<List<Integer>> ans, List<Integer> temp, boolean[] freq) {
            if (temp.size() == nums.length) {
                ans.add(new ArrayList<>(temp));
                return;
            }

            for (int i = 0; i < nums.length; i++) {
                if (!freq[i]) // if it has not been selected
                {
                    freq[i] = true;
                    temp.add(nums[i]); // ad it to the DS
                    backtrack(nums, ans, temp, freq); // no change in params
                    temp.remove(temp.size() - 1);// remove the last element
                    // make it false, since it is not picked now and can be re picked again
                    freq[i] = false;
                }
            }
        }
    }

    // ^ N queens Problem

    // ^ Sudoku Solver

    // ^ M coloring Problem

    // ^ Rat in a Maze

    // ^ Word Break (print all ways)
}
