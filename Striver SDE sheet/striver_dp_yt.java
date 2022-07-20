import java.util.*;

class striver_dp_yt {
    // ^ Frog jump ( lecture 3 )

    // There is a frog on the 1st step of an N stairs long staircase. The frog wants
    // to reach the Nth stair. HEIGHT[i] is the height of the (i+1)th stair.If Frog
    // jumps from ith to jth stair, the energy lost in the jump is given by
    // |HEIGHT[i-1] - HEIGHT[j-1] |.In the Frog is on ith staircase, he can jump
    // either to (i+1)th stair or to (i+2)th stair. Your task is to find the minimum
    // total energy used by the frog to reach from 1st stair to Nth stair.

    // For Example
    // If the given ‘HEIGHT’ array is [10,20,30,10], the answer 20 as the frog can
    // jump from 1st stair to 2nd stair (|20-10| = 10 energy lost) and then a jump
    // from 2nd stair to last stair (|10-20| = 10 energy lost). So, the total energy
    // lost is 20.

    // ~ naive solution using only recursion
    public class frog_jump_recursion {
        public int frogJump(int n, int heights[]) {
            return fun(n - 1, heights);
        }

        // trying out all possible ways and taking the best
        private int fun(int n, int[] heights) {
            if (n == 0)
                return 0;
            int left = fun(n - 1, heights) + Math.abs(heights[n] - heights[n - 1]);
            int right = Integer.MAX_VALUE;
            if (n > 1) {
                right = fun(n - 2, heights) + Math.abs(heights[n] - heights[n - 2]);
            }
            return Math.min(left, right);
        }
    }

    // ~ using memoization technique so already computed values are not computed
    // again

    public class frog_jump_memoization {
        public int frogJump(int n, int heights[]) {
            int[] dp = new int[n + 1];
            Arrays.fill(dp, -1);
            return fun(n - 1, heights, dp);
        }

        private int fun(int n, int[] heights, int[] dp) {
            if (n == 0)
                return 0;
            if (dp[n] != -1)
                return dp[n];
            int left = fun(n - 1, heights, dp) + Math.abs(heights[n] - heights[n - 1]);
            int right = Integer.MAX_VALUE;
            if (n > 1) {
                right = fun(n - 2, heights, dp) + Math.abs(heights[n] - heights[n - 2]);
            }
            return dp[n] = Math.min(left, right);
        }
    }

    // ~ using tabulation, just replace the function call with dp array
    // tabulation uses bottom up approach
    public class frog_jump_tabulation {
        public int frogJump(int n, int heights[]) {
            int[] dp = new int[n + 1];
            Arrays.fill(dp, 0);
            for (int i = 1; i < n; i++) {
                int left = dp[i - 1] + Math.abs(heights[i] - heights[i - 1]);
                int right = Integer.MAX_VALUE;
                if (i > 1) {
                    right = dp[i - 2] + Math.abs(heights[i] - heights[i - 2]);
                }
                dp[i] = Math.min(left, right);
            }
            return dp[n - 1];
        }
    }

    // instead of using a dp array, use 2 variables prev nd prev2
    public class frog_jump_space_optimized {
        public int frogJump(int n, int heights[]) {
            int prev = 0, prev2 = 0;
            for (int i = 1; i < n; i++) {
                int left = prev + Math.abs(heights[i] - heights[i - 1]);
                int right = Integer.MAX_VALUE;
                if (i > 1) {
                    right = prev2 + Math.abs(heights[i] - heights[i - 2]);
                }
                int curr = Math.min(left, right);
                prev2 = prev;
                prev = curr;
            }
            return prev;
        }
    }

    // ^ Max sum of non adjacent elements

    // You are given an array/list of ‘N’ integers. You are supposed to return the
    // maximum sum of the subsequence with the constraint that no two elements are
    // adjacent in the given array/list.

    // Sample Input 1:

    // 1 2 4

    // 2 1 4 9

    // Sample Output 1:
    // 5
    // 11

    // ~ Basic Recursion approach

    public class max_sum_non_adjacent {
        public int maximumNonAdjacentSum(ArrayList<Integer> nums) {
            int n = nums.size();
            return f(n - 1, nums);
        }

        public int f(int index, ArrayList<Integer> nums) {
            if (index == 0)
                return nums.get(index);
            if (index < 0)
                return 0;

            int pick = nums.get(index) + f(index - 2, nums);
            int not_pick = 0 + f(index - 1, nums);

            return Math.max(pick, not_pick);
        }
    }

    // ~ Applying Memoization

    public class max_sum_non_adjacent_memoization {
        public int maximumNonAdjacentSum(ArrayList<Integer> nums) {
            int n = nums.size();
            int[] dp = new int[n];
            Arrays.fill(dp, -1);
            return f(n - 1, nums, dp);
        }

        public int f(int index, ArrayList<Integer> nums, int[] dp) {
            if (index == 0)
                return nums.get(index);
            if (index < 0)
                return 0;

            // skipping overlapping sub-problems
            if (dp[index] != -1)
                return dp[index];
            int pick = nums.get(index) + f(index - 2, nums, dp);
            int not_pick = 0 + f(index - 1, nums, dp);

            return dp[index] = Math.max(pick, not_pick);
        }
    }

    // ~ Tabulation with space optimization

    public class max_sum_non_adjacent_tabulation {
        public int maximumNonAdjacentSum(ArrayList<Integer> nums) {
            int n = nums.size();
            int prev = nums.get(0);
            int prev2 = 0;

            // replace the dp array
            for (int i = 1; i < n; i++) {
                int take = nums.get(i);
                if (i > 1) {
                    take += prev2;
                }
                int not_take = 0 + prev;

                int curri = Math.max(take, not_take);

                prev2 = prev;
                prev = curri;
            }
            return prev; // because we would have returned dp[n-1]
        }
    }

    // ^ House Robber
    // Same as above
    // Refer Notebook

    // ^ Ninja's Training

    // Ninja is planing this ‘N’ days-long training schedule. Each day, he can
    // perform any one of these three activities. (Running, Fighting Practice or
    // Learning New Moves). Each activity has some merit points on each day. As
    // Ninja has to improve all his skills, he can’t do the same activity in two
    // consecutive days. Can you help Ninja find out the maximum merit points Ninja
    // can earn?

    // You are given a 2D array of size N*3 ‘POINTS’ with the points corresponding
    // to each day and activity. Your task is to calculate the maximum number of
    // merit points that Ninja can earn.

    // ~ Basic Recursion Approach

    public class Ninja_training {

        public int f(int day, int last, int[][] points) {
            if (day == 0) {
                int maxi = 0;
                for (int i = 0; i < 3; i++) {
                    if (i != last)
                        maxi = Math.max(maxi, points[0][i]);
                }
                return maxi;
            }

            int maxi = 0;
            for (int i = 0; i < 3; i++) {
                if (i != last) {
                    int point = points[day][i] + f(day - 1, i, points);
                    maxi = Math.max(maxi, point);
                }
            }
            return maxi;
        }

        public int ninjaTraining(int n, int points[][]) {

            return f(n - 1, 3, points);
        }
    }

    // ~ Memoization

    public class Ninja_training_memoization {

        public int f(int day, int last, int[][] points, int[][] dp) {

            if (dp[day][last] != -1)
                return dp[day][last];

            if (day == 0) {
                int maxi = 0;
                for (int i = 0; i < 3; i++) {
                    if (i != last)
                        maxi = Math.max(maxi, points[0][i]);
                }
                return dp[day][last] = maxi;
            }

            int maxi = 0;
            for (int i = 0; i < 3; i++) {
                if (i != last) {
                    int point = points[day][i] + f(day - 1, i, points, dp);
                    maxi = Math.max(maxi, point);
                }
            }
            return dp[day][last] = maxi;
        }

        public int ninjaTraining(int n, int points[][]) {

            int[][] dp = new int[n][4];
            for (int[] row : dp) {
                Arrays.fill(row, -1);
            }
            return f(n - 1, 3, points, dp);
        }
    }

    // ~ Tabulation

    public class Ninja_training_tabulation {

        public int ninjaTraining(int n, int points[][]) {

            int[][] dp = new int[n][4];
            for (int[] row : dp) {
                Arrays.fill(row, -1);
            }

            dp[0][0] = Math.max(points[0][1], points[0][2]);
            dp[0][1] = Math.max(points[0][0], points[0][2]);
            dp[0][2] = Math.max(points[0][0], points[0][0]);
            dp[0][3] = Math.max(points[0][1], Math.max(points[0][2], points[0][2]));

            for (int day = 1; day < n; day++) {
                for (int last = 0; last < 4; last++) {
                    // copy paste the recursion from here
                    dp[day][last] = 0;

                    for (int task = 0; task < 3; task++) {
                        if (task != last) {
                            int point = points[day][task] + dp[day - 1][task];
                            dp[day][last] = Math.max(dp[day][last], point);
                        }
                    }
                }
            }
            return dp[n - 1][3]; // just as in recursion we were calling the recursive function with the
                                 // arguments
        }
    }

    // ~ Space optimized Tabulation

    public class Ninja_training_space_optimized {

        public int ninjaTraining(int n, int points[][]) {

            int[] dp = new int[4];
            Arrays.fill(dp, 0);

            dp[0] = Math.max(points[0][1], points[0][2]);
            dp[1] = Math.max(points[0][0], points[0][2]);
            dp[2] = Math.max(points[0][0], points[0][0]);
            dp[3] = Math.max(points[0][1], Math.max(points[0][2], points[0][2]));

            for (int day = 1; day < n; day++) {
                int[] temp = new int[4];
                Arrays.fill(temp, 0);
                for (int last = 0; last < 4; last++) {
                    // copy paste the recursion from here
                    temp[last] = 0;

                    for (int task = 0; task < 3; task++) {
                        if (task != last) {
                            int point = points[day][task] + dp[task];
                            temp[last] = Math.max(temp[last], point);
                        }
                    }
                }
                dp = temp;
            }
            return dp[3]; // just as in recursion we were calling the recursive function with there
                          // arguments
        }
    }

    // ^ Grid unique Paths

    // You are present at point ‘A’ which is the top-left cell of an M X N matrix,
    // your destination is point ‘B’, which is the bottom-right cell of the same
    // matrix. Your task is to find the total number of unique paths from point ‘A’
    // to point ‘B’.In other words, you will be given the dimensions of the matrix
    // as integers ‘M’ and ‘N’, your task is to find the total number of unique
    // paths from the cell MATRIX[0][0] to MATRIX['M' - 1]['N' - 1].

    // To traverse in the matrix, you can either move Right or Down at each step.
    // For example in a given point MATRIX[i] [j], you can move to either MATRIX[i +
    // 1][j] or MATRIX[i][j + 1].

    // ~ Basic Recursive Approach

    public class grid_unique_paths_recursion {
        public int uniquePaths(int m, int n) {
            return f(m - 1, n - 1);
        }

        public int f(int i, int j) {
            if (i == 0 && j == 0)
                return 1;
            if (i < 0 || j < 0)
                return 0;

            int up = f(i - 1, j);
            int down = f(i, j - 1);

            return up + down;
        }
    }

    // ~ Memoization

    public class grid_unique_paths_memoization {
        public int uniquePaths(int m, int n) {
            int[][] dp = new int[m][n];
            for (int[] row : dp) {
                Arrays.fill(row, -1);
            }
            return f(m - 1, n - 1, dp);
        }

        public int f(int i, int j, int[][] dp) {
            if (i == 0 && j == 0)
                return 1;
            if (i < 0 || j < 0)
                return 0;

            if (dp[i][j] != -1)
                return dp[i][j];

            int up = f(i - 1, j, dp);
            int down = f(i, j - 1, dp);

            return dp[i][j] = up + down;
        }
    }

    // ~ Tabulation

    public class grid_unique_paths_tabulation {
        public int uniquePaths(int m, int n) {
            int[][] dp = new int[m][n];

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == 0 && j == 0)
                        dp[i][j] = 1;
                    else {
                        int up = 0, left = 0;
                        if (i > 0)
                            up = dp[i - 1][j];
                        if (j > 0)
                            left = dp[i][j - 1];

                        dp[i][j] = up + left;
                    }
                }
            }
            return dp[m - 1][n - 1];
        }
    }

    // ~ Space optimized Tabulation

    public class grid_unique_paths_space_optimized {
        public int uniquePaths(int m, int n) {

            int[] prev = new int[n]; // dp is prev row

            for (int i = 0; i < m; i++) {
                int[] curr = new int[n]; // curr is current row
                for (int j = 0; j < n; j++) {
                    if (i == 0 && j == 0)
                        curr[j] = 1;
                    else {
                        // here dp[i] is current and dp[i-1] is prev
                        int up = 0, left = 0;
                        if (i > 0)
                            up = prev[j];
                        if (j > 0)
                            left = curr[j - 1];

                        curr[j] = up + left;
                    }
                }
                prev = curr;
            }
            return prev[n - 1];
        }
    }

    // ^ Grid unique paths 2
    // same previous problem with some obstacles in between

    // ^ Minimum Path Sum
    // Ninja land is a country in the shape of a 2-Dimensional grid 'GRID', with 'N'
    // rows and 'M' columns. Each point in the grid has some cost associated with
    // it.
    // Find a path from top left i.e. (0, 0) to the bottom right i.e. ('N' - 1, 'M'
    // - 1) which minimizes the sum of the cost of all the numbers along the path.
    // You need to tell the minimum sum of that path.

    // ~ Basic Recursion Approach

    public class min_path_sum {
        public int minSumPath(int[][] grid) {
            int m = grid.length;
            int n = grid[0].length;
            return f(m - 1, n - 1, grid);
        }

        public int f(int i, int j, int[][] grid) {
            if (i == 0 && j == 0)
                return grid[i][j];
            if (i < 0 || j < 0)
                return (int) Math.pow(10, 9);

            int up = grid[i][j] + f(i - 1, j, grid);
            int left = grid[i][j] + f(i, j - 1, grid);

            return Math.min(up, left);
        }
    }

    // ~ Memoization

    public class min_path_sum_memoization {
        public int minSumPath(int[][] grid) {
            int m = grid.length;
            int n = grid[0].length;
            int[][] dp = new int[m][n];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(m - 1, n - 1, grid, dp);
        }

        public int f(int i, int j, int[][] grid, int[][] dp) {
            if (i == 0 && j == 0)
                return grid[i][j];
            if (i < 0 || j < 0)
                return (int) Math.pow(10, 9);

            if (dp[i][j] != -1)
                return dp[i][j];

            int up = grid[i][j] + f(i - 1, j, grid, dp);
            int left = grid[i][j] + f(i, j - 1, grid, dp);

            return dp[i][j] = Math.min(up, left);
        }
    }

    // ~ Tabulation

    public class min_path_sum_tabulation {
        public int minSumPath(int[][] grid) {
            int m = grid.length;
            int n = grid[0].length;
            int[][] dp = new int[m][n];
            for (int[] row : dp)
                Arrays.fill(row, 0);

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == 0 && j == 0)
                        dp[i][j] = grid[i][j];
                    else {
                        int up = grid[i][j], left = grid[i][j];
                        if (i > 0)
                            up += dp[i - 1][j];
                        else
                            up = Integer.MAX_VALUE;

                        if (j > 0)
                            left += dp[i][j - 1];
                        else
                            left = Integer.MAX_VALUE;

                        dp[i][j] = Math.min(up, left);
                    }
                }
            }
            return dp[m - 1][n - 1];
        }
    }

    // ~ Space optimized

    public class min_path_sum_space_optimized {
        public int minSumPath(int[][] grid) {
            int m = grid.length;
            int n = grid[0].length;
            int[] prev = new int[n];

            // in tabulation code dp[i][j] is curr[j] and dp[i-1][j] is prev[j]
            // we replace matrix with prev array

            for (int i = 0; i < m; i++) {
                int[] curr = new int[n];
                for (int j = 0; j < n; j++) {
                    if (i == 0 && j == 0)
                        curr[j] = grid[i][j];
                    else {
                        int up = grid[i][j], left = grid[i][j];
                        // requiring previous row's j column
                        if (i > 0)
                            up += prev[j];
                        else
                            up = Integer.MAX_VALUE;

                        // requiring current row's j-1 column
                        if (j > 0)
                            left += curr[j - 1];
                        else
                            left = Integer.MAX_VALUE;

                        curr[j] = Math.min(up, left);
                    }
                }
                prev = curr;
            }
            return prev[n - 1];
        }
    }
}