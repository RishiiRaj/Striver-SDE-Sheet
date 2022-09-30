import java.util.*;

class striver_dp_yt {

    // ! lecture 3

    // ^ Frog jump

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

    // ~ Recursion

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

    // ~ Memoization

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

    // ~ Tabulation

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

    // ~ Space Optimized

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

    // ! Lecture 5

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

    // ! Lecture 6

    // ^ House Robber
    // Same as above
    // Refer Notebook

    // ! Lecture 7

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

    // ! Lecture 8

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

    // ! Lecture 9

    // ^ Grid unique paths 2
    // same previous problem with some obstacles in between

    // ! Lecture 10

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

    // ! Lecture 11

    // ^ Triangle

    // You are given a triangular array/list 'TRIANGLE'. Your task is to return the
    // minimum path sum to reach from the top to the bottom row.
    // The triangle array will have N rows and the i-th row, where 0 <= i < N will
    // have i + 1 elements.
    // You can move only to the adjacent number of row below each step. For example,
    // if you are at index j in row i, then you can move to i or i + 1 index in row
    // j + 1 in each step.
    // For Example :
    // If the array given is 'TRIANGLE' = [[1], [2,3], [3,6,7], [8,9,6,1]] the
    // triangle array will look like:

    // 1
    // 2,3
    // 3,6,7
    // 8,9,6,10

    // For the given triangle array the minimum sum path would be 1->2->3->8. Hence
    // the answer would be 14.

    // ~ Recursion

    public class triangle_recursion {
        public int minimumPathSum(int[][] triangle, int n) {
            return f(0, 0, triangle, n);
        }

        public int f(int i, int j, int[][] triangle, int n) {
            if (i == n - 1) {
                return triangle[i][j];
            }

            int down = triangle[i][j] + f(i + 1, j, triangle, n);
            int diag = triangle[i][j] + f(i + 1, j + 1, triangle, n);

            return Math.min(down, diag);
        }
    }

    // ~ Memoization

    public class triangle_memoization {
        public int minimumPathSum(int[][] triangle, int n) {
            int[][] dp = new int[n][n];
            for (int[] row : dp)
                Arrays.fill(row, -1);

            return f(0, 0, triangle, n, dp);
        }

        public int f(int i, int j, int[][] triangle, int n, int[][] dp) {
            if (i == n - 1) {
                return triangle[i][j];
            }

            if (dp[i][j] != -1)
                return dp[i][j];

            int down = triangle[i][j] + f(i + 1, j, triangle, n, dp);
            int diag = triangle[i][j] + f(i + 1, j + 1, triangle, n, dp);

            return dp[i][j] = Math.min(down, diag);
        }
    }

    // ~ Tabulation

    public class triangle_tabulation {
        public int minimumPathSum(int[][] triangle, int n) {
            int[][] dp = new int[n][n];
            for (int[] row : dp)
                Arrays.fill(row, 0);

            for (int j = 0; j < n; j++)
                dp[n - 1][j] = triangle[n - 1][j];

            for (int i = n - 2; i >= 0; i--) {
                for (int j = i; j >= 0; j--) {
                    int down = triangle[i][j] + dp[i + 1][j];
                    int diag = triangle[i][j] + dp[i + 1][j + 1];

                    dp[i][j] = Math.min(down, diag);
                }
            }
            return dp[0][0];
        }
    }

    // ~ Space Optimization

    public class triangle_space_optimized {
        public int minimumPathSum(int[][] triangle, int n) {
            int[] front = new int[n];
            int[] curr = new int[n];

            for (int j = 0; j < n; j++)
                front[j] = triangle[n - 1][j];

            for (int i = n - 2; i >= 0; i--) {
                for (int j = i; j >= 0; j--) {
                    int down = triangle[i][j] + front[j];
                    int diag = triangle[i][j] + front[j + 1];
                    curr[j] = Math.min(down, diag);
                }
                front = curr;
            }
            return front[0];
        }
    }

    // ! Lecture 12

    // ^ Maximum Path Sum

    // We are given an ‘N*M’ matrix. We need to find the maximum path sum from any
    // cell of the first row to any cell of the last row.

    // At every cell we can move in three directions: to the bottom cell (↓), to the
    // bottom-right cell(↘), or to the bottom-left cell(↙).

    // ~ Recursion

    public class max_path_sum_recursion {
        public int getMaxPathSum(int[][] matrix) {
            int n = matrix.length;
            int m = matrix[0].length;

            // finding maximum of that row
            int maxi = Integer.MIN_VALUE;
            for (int j = 0; j < m; j++) {
                int ans = f(n - 1, j, m, matrix);
                maxi = Math.max(maxi, ans);
            }
            return maxi;
        }

        public int f(int i, int j, int m, int[][] a) {
            if (j < 0 || j >= m)
                return (int) Math.pow(-10, 9); // a very small value to omit the Max func
            if (i == 0)
                return a[0][j];

            int straight = a[i][j] + f(i - 1, j, m, a);
            int left_diag = a[i][j] + f(i - 1, j - 1, m, a);
            int right_diag = a[i][j] + f(i - 1, j + 1, m, a);

            return Math.max(straight, Math.max(left_diag, right_diag));
        }
    }

    // ~ Memoization

    public class max_path_sum_memoization {
        public int getMaxPathSum(int[][] matrix) {
            int n = matrix.length;
            int m = matrix[0].length;

            int[][] dp = new int[n][m];
            for (int[] row : dp)
                Arrays.fill(row, -1);

            int maxi = Integer.MIN_VALUE;
            for (int j = 0; j < m; j++) {
                int ans = f(n - 1, j, m, matrix, dp);
                maxi = Math.max(maxi, ans);
            }
            return maxi;
        }

        public int f(int i, int j, int m, int[][] a, int[][] dp) {
            if (j < 0 || j >= m)
                return (int) Math.pow(-10, 9);
            if (i == 0)
                return a[0][j];

            if (dp[i][j] != -1)
                return dp[i][j];

            int straight = a[i][j] + f(i - 1, j, m, a, dp);
            int left_diag = a[i][j] + f(i - 1, j - 1, m, a, dp);
            int right_diag = a[i][j] + f(i - 1, j + 1, m, a, dp);

            return dp[i][j] = Math.max(straight, Math.max(left_diag, right_diag));
        }
    }

    // ~ tabulation

    public class max_path_sum_tabulation {
        public int getMaxPathSum(int[][] matrix) {
            int n = matrix.length;
            int m = matrix[0].length;
            int[][] dp = new int[n][m];

            for (int j = 0; j < m; j++)
                dp[0][j] = matrix[0][j];

            for (int i = 1; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    int straight = matrix[i][j] + dp[i - 1][j];

                    int ld = matrix[i][j];
                    if (j - 1 >= 0)
                        ld = matrix[i][j] + dp[i - 1][j - 1];
                    else
                        ld = (int) Math.pow(-10, 9);

                    int rd = matrix[i][j];
                    if (j + 1 < m)
                        rd = matrix[i][j] + dp[i - 1][j + 1];
                    else
                        rd = (int) Math.pow(-10, 9);

                    dp[i][j] = Math.max(straight, Math.max(ld, rd));
                }
            }

            int maxi = Integer.MIN_VALUE;
            for (int j = 0; j < m; j++) {
                int ans = dp[n - 1][j];// f(n-1,j,m,matrix,dp);
                maxi = Math.max(maxi, ans);
            }
            return maxi;
        }
    }

    // ~ Space optimization

    public class max_path_sum_space_optimized {
        public int getMaxPathSum(int[][] matrix) {
            int n = matrix.length;
            int m = matrix[0].length;
            // int[][] dp = new int[n][m];
            int[] prev = new int[m];
            int[] curr = new int[m];

            for (int j = 0; j < m; j++)
                prev[j] = matrix[0][j];

            for (int i = 1; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    int straight = matrix[i][j] + prev[j];

                    int ld = matrix[i][j];
                    if (j - 1 >= 0)
                        ld = matrix[i][j] + prev[j - 1];
                    else
                        ld = (int) Math.pow(-10, 9);

                    int rd = matrix[i][j];
                    if (j + 1 < m)
                        rd = matrix[i][j] + prev[j + 1];
                    else
                        rd = (int) Math.pow(-10, 9);

                    curr[j] = Math.max(straight, Math.max(ld, rd));
                }
                prev = curr;
            }

            int maxi = Integer.MIN_VALUE;
            for (int j = 0; j < m; j++) {
                int ans = prev[j];// f(n-1,j,m,matrix,dp);
                maxi = Math.max(maxi, ans);
            }
            return maxi;
        }
    }

    // ! Lecture 13

    // ^ Ninja and his Friends

    // We are given an ‘N*M’ matrix. Every cell of the matrix has some chocolates on
    // it, mat[i][j] gives us the number of chocolates. We have two friends ‘Alice’
    // and ‘Bob’. initially, Alice is standing on the cell(0,0) and Bob is standing
    // on the cell(0, M-1). Both of them can move only to the cells below them in
    // these three directions: to the bottom cell (↓), to the bottom-right cell(↘),
    // or to the bottom-left cell(↙).

    // When Alice and Bob visit a cell, they take all the chocolates from that cell
    // with them. It can happen that they visit the same cell, in that case, the
    // chocolates need to be considered only once.

    // They cannot go out of the boundary of the given matrix, we need to return the
    // maximum number of chocolates that Bob and Alice can together collect.

    // ~ Recursion

    public class ninja_and_friends_recursion {
        public int maximumChocolates(int r, int c, int[][] grid) {
            return f(0, 0, c - 1, grid);
        }

        public int f(int i, int j1, int j2, int[][] a) {
            int m = a[0].length;
            int n = a.length;
            // out of bounds
            if (j1 < 0 || j1 >= m || j2 < 0 || j2 >= m)
                return (int) Math.pow(-10, 9);

            // reached the last row
            if (i == n - 1) {
                // if both end up at the same cell
                if (j1 == j2)
                    return a[i][j1]; // return any one
                else
                    return a[i][j1] + a[i][j2];
            }
            // explore all paths of alice and bob simultaneously
            int maxi = Integer.MIN_VALUE;

            for (int dj1 = -1; dj1 <= 1; dj1++) {
                for (int dj2 = -1; dj2 <= 1; dj2++) {
                    // both at same cell
                    int value = 0;
                    if (j1 == j2)
                        value = a[i][j1];
                    else
                        value = a[i][j1] + a[i][j2];
                    value += f(i + 1, j1 + dj1, j2 + dj2, a);
                    maxi = Math.max(maxi, value);
                }
            }
            return maxi;
        }
    }

    // ~ Memoization

    public class ninja_and_friends_memoization {
        public int maximumChocolates(int r, int c, int[][] grid) {
            int[][][] dp = new int[r][c][c];
            for (int row1[][] : dp)
                for (int row2[] : row1)
                    Arrays.fill(row2, -1);

            return f(0, 0, c - 1, grid, dp);
        }

        public int f(int i, int j1, int j2, int[][] a, int[][][] dp) {
            int m = a[0].length;
            int n = a.length;
            // out of bounds
            if (j1 < 0 || j1 >= m || j2 < 0 || j2 >= m)
                return (int) Math.pow(-10, 9);

            // reached the last row
            if (i == n - 1) {
                // if both end up at the same cell
                if (j1 == j2)
                    return a[i][j1]; // return any one
                else
                    return a[i][j1] + a[i][j2];
            }
            if (dp[i][j1][j2] != -1)
                return dp[i][j1][j2];
            // explore all paths of alice and bob simultaneously
            int maxi = Integer.MIN_VALUE;

            for (int dj1 = -1; dj1 <= 1; dj1++) {
                for (int dj2 = -1; dj2 <= 1; dj2++) {
                    // both at same cell
                    int value = 0;
                    if (j1 == j2)
                        value = a[i][j1];
                    else
                        value = a[i][j1] + a[i][j2];
                    value += f(i + 1, j1 + dj1, j2 + dj2, a, dp);
                    maxi = Math.max(maxi, value);
                }
            }
            return dp[i][j1][j2] = maxi;
        }
    }

    // ~ Tabulation

    public class ninja_and_friends_tabulation {
        public int maximumChocolates(int n, int m, int[][] grid) {
            int[][][] dp = new int[n][m][m];

            // Base Case
            for (int j1 = 0; j1 < m; j1++) {
                for (int j2 = 0; j2 < m; j2++) {
                    if (j1 == j2)
                        dp[n - 1][j1][j2] = grid[n - 1][j1];
                    else
                        dp[n - 1][j1][j2] = grid[n - 1][j1] + grid[n - 1][j2];
                }
            }

            for (int i = n - 2; i >= 0; i--) {
                for (int j1 = 0; j1 < m; j1++) {
                    for (int j2 = 0; j2 < m; j2++) {
                        // just copy the recurrence
                        int maxi = Integer.MIN_VALUE;

                        for (int dj1 = -1; dj1 <= 1; dj1++) {
                            for (int dj2 = -1; dj2 <= 1; dj2++) {
                                // both at same cell
                                int value = 0;
                                if (j1 == j2)
                                    value = grid[i][j1];
                                else
                                    value = grid[i][j1] + grid[i][j2];

                                // checking out of bounds condition
                                if (j1 + dj1 >= 0 && j1 + dj1 < m && j2 + dj2 >= 0 && j2 + dj2 < m)
                                    value += dp[i + 1][j1 + dj1][j2 + dj2];
                                maxi = Math.max(maxi, value);
                            }
                        }
                        dp[i][j1][j2] = maxi;
                    }
                }
            }
            return dp[0][0][m - 1];
        }
    }

    // ~ Space optimization

    public class ninja_and_friends_space_optimized {
        public int maximumChocolates(int n, int m, int[][] grid) {

            int[][] front = new int[m][m];
            int[][] cur = new int[m][m];

            // change dp[n-1] to front
            for (int j1 = 0; j1 < m; j1++) {
                for (int j2 = 0; j2 < m; j2++) {
                    if (j1 == j2)
                        front[j1][j2] = grid[n - 1][j1];
                    else
                        front[j1][j2] = grid[n - 1][j1] + grid[n - 1][j2];
                }
            }

            // Outer Nested Loops for traversing DP Array
            for (int i = n - 2; i >= 0; i--) {
                for (int j1 = 0; j1 < m; j1++) {
                    for (int j2 = 0; j2 < m; j2++) {

                        int maxi = Integer.MIN_VALUE;

                        // change dp[i+1] to front
                        // Inner nested loops to try out 9 options
                        for (int dj1 = -1; dj1 <= 1; dj1++) {
                            for (int dj2 = -1; dj2 <= 1; dj2++) {

                                int ans;

                                if (j1 == j2)
                                    ans = grid[i][j1];
                                else
                                    ans = grid[i][j1] + grid[i][j2];

                                if ((j1 + dj1 < 0 || j1 + dj1 >= m) ||
                                        (j2 + dj2 < 0 || j2 + dj2 >= m))

                                    ans += (int) Math.pow(-10, 9);
                                else
                                    ans += front[j1 + dj1][j2 + dj2];

                                maxi = Math.max(ans, maxi);
                            }
                        }
                        // change dp[0] to curr
                        cur[j1][j2] = maxi;
                    }
                }

                for (int a = 0; a < m; a++) {
                    front[a] = (int[]) (cur[a].clone());
                }
            }

            return front[0][m - 1];
        }
    }

    // ! Lecture 14

    // ^ Subset/Subsequence sum equals target

    // We are given an array ‘ARR’ with N positive integers. We need to find if
    // there is a subset in “ARR” with a sum equal to K. If there is, return true
    // else return false.

    // ~ Recursion

    public class subset_sum_equals_K_recursion {
        public boolean subsetSumToK(int n, int k, int arr[]) {
            return f(n - 1, k, arr);
        }

        public boolean f(int index, int target, int[] a) {
            if (target == 0)
                return true;
            if (index == 0)
                return (a[0] == target);

            boolean not_take = f(index - 1, target, a);
            boolean take = false;
            if (target >= a[index])
                take = f(index - 1, target - a[index], a);

            return (take | not_take);
        }
    }

    // ~ Memoization

    public class subset_sum_equals_K_memoization {
        public boolean subsetSumToK(int n, int k, int arr[]) {
            int[][] dp = new int[n + 1][k + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, k, arr, dp);
        }

        public boolean f(int index, int target, int[] a, int[][] dp) {
            if (target == 0)
                return true;
            if (index == 0)
                return (a[0] == target);

            if (dp[index][target] != -1)
                return dp[index][target] == 0 ? false : true;
            boolean not_take = f(index - 1, target, a, dp);
            boolean take = false;
            if (target >= a[index])
                take = f(index - 1, target - a[index], a, dp);

            dp[index][target] = not_take || take ? 1 : 0;
            return (take | not_take);
        }
    }

    // ~ Tabulation

    public class subset_sum_equals_K_tabulation {
        public boolean subsetSumToK(int n, int k, int arr[]) {
            boolean[][] dp = new boolean[n + 1][k + 1];
            // first base case
            for (int i = 0; i < n; i++) {
                dp[i][0] = true;
            }
            // second base case
            if (arr[0] < k)
                dp[0][arr[0]] = true;

            for (int index = 1; index < n; index++) {
                for (int target = 1; target <= k; target++) {
                    // just copy the recurrence and change recursion call to dp matrix
                    boolean not_take = dp[index - 1][target];
                    boolean take = false;
                    if (target >= arr[index])
                        take = dp[index - 1][target - arr[index]];

                    dp[index][target] = (take | not_take);
                }
            }
            // how did we call the f function
            return dp[n - 1][k];
        }
    }

    // ~ Space Optimization

    public class subset_sum_equals_K_space_optimized {
        public boolean subsetSumToK(int n, int k, int arr[]) {
            boolean[] prev = new boolean[k + 1];

            // second base case
            if (arr[0] <= k)
                prev[arr[0]] = true;

            for (int index = 1; index < n; index++) {
                boolean[] curr = new boolean[k + 1];
                curr[0] = true;
                for (int target = 1; target <= k; target++) {
                    // just copy the recurrence and change recursion call to dp matrix
                    boolean not_take = prev[target];
                    boolean take = false;
                    if (target >= arr[index])
                        take = prev[target - arr[index]];

                    curr[target] = (take || not_take);
                }
                prev = curr;
            }
            // how did we call the f function
            return prev[k];
        }
    }

    // ! Lecture 15

    // ^ Partition Equal Subset Sum

    // We are given an array ‘ARR’ with N positive integers. We need to find if we
    // can partition the array into two subsets such that the sum of elements of
    // each subset is equal to the other.

    public class partition_equak_subset_sum {
        public boolean canPartition(int[] arr, int n) {
            int sum = 0;
            for (int i : arr)
                sum += i;
            if (sum % 2 != 0)
                return false;
            else
                return subsetSumToK(n, sum / 2, arr);
        }

        // previous lecture function
        public boolean subsetSumToK(int n, int k, int arr[]) {
            boolean[] prev = new boolean[k + 1];

            // second base case
            if (arr[0] <= k)
                prev[arr[0]] = true;

            for (int index = 1; index < n; index++) {
                boolean[] curr = new boolean[k + 1];
                curr[0] = true;
                for (int target = 1; target <= k; target++) {
                    // just copy the recurrence and change recursion call to dp matrix
                    boolean not_take = prev[target];
                    boolean take = false;
                    if (target >= arr[index])
                        take = prev[target - arr[index]];

                    curr[target] = (take || not_take);
                }
                prev = curr;
            }
            // how did we call the f function
            return prev[k];
        }
    }

    // ! Lecture 16

    // ^ Partition a set into two subsets with minimum absolute difference

    public class min_abs_subset_sum_difference {
        public int minSubsetSumDifference(int[] arr, int n) {
            int sum = 0;
            for (int i : arr)
                sum += i;
            int k = sum;
            // * just using the L-14 tabulation function
            boolean[][] dp = new boolean[n + 1][k + 1];
            // first base case
            for (int i = 0; i < n; i++) {
                dp[i][0] = true;
            }
            // second base case
            if (arr[0] <= k)
                dp[0][arr[0]] = true;

            for (int index = 1; index < n; index++) {
                for (int target = 1; target <= k; target++) {
                    // just copy the recurrence and change recursion call to dp matrix
                    boolean not_take = dp[index - 1][target];
                    boolean take = false;
                    if (target >= arr[index])
                        take = dp[index - 1][target - arr[index]];

                    dp[index][target] = (take | not_take);
                }
            }

            // * The Extra part in this question
            // dp[n-1][col-> 0......total_sum] will give the ans
            int mini = Integer.MAX_VALUE;
            for (int s1 = 0; s1 <= sum / 2; s1++) {
                if (dp[n - 1][s1] == true) {
                    int s2 = sum - s1;
                    mini = Math.min(mini, Math.abs(s2 - s1));
                }
            }
            return mini;
        }
    }

    // ! Lecture 17

    // ^ Count subsets with sum K

    // ~ Recursion

    public class count_subsets_with_sum_K_recursion {
        public int findWays(int num[], int tar) {
            int n = num.length;
            return f(n - 1, tar, num);
        }

        public int f(int index, int sum, int[] a) {
            if (sum == 0)
                return 1;
            if (index == 0) {
                if (a[index] == sum)
                    return 1;
                else
                    return 0;
            }
            int not_pick = f(index - 1, sum, a);
            int pick = 0;
            if (sum >= a[index])
                pick = f(index - 1, sum - a[index], a);
            return pick + not_pick;
        }
    }

    // ~ Memoization

    public class count_subsets_with_sum_K_memoization {
        public int findWays(int num[], int tar) {
            int n = num.length;
            int[][] dp = new int[n][tar + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, tar, num, dp);
        }

        public int f(int index, int sum, int[] a, int[][] dp) {
            if (sum == 0)
                return 1;
            if (index == 0) {
                if (a[index] == sum)
                    return 1;
                else
                    return 0;
            }
            if (dp[index][sum] != -1)
                return dp[index][sum];

            int not_pick = f(index - 1, sum, a, dp);
            int pick = 0;
            if (sum >= a[index])
                pick = f(index - 1, sum - a[index], a, dp);
            return dp[index][sum] = pick + not_pick;
        }
    }

    // ~ Tabulation

    public class count_subsets_with_sum_K_tabulation {
        public int findWays(int num[], int tar) {
            int n = num.length;
            int[][] dp = new int[n][tar + 1];

            // base cases
            for (int i = 0; i < n; i++) {
                dp[i][0] = 1;
            }
            if (tar >= num[0])
                dp[0][num[0]] = 1;

            // look at changing parameters and create nested loops
            for (int index = 1; index < n; index++) {
                for (int sum = 1; sum <= tar; sum++) {
                    // copy the recurrence
                    int not_pick = dp[index - 1][sum];
                    int pick = 0;
                    if (sum >= num[index])
                        pick = dp[index - 1][sum - num[index]];
                    dp[index][sum] = pick + not_pick;
                }
            }
            return dp[n - 1][tar];
        }
    }

    // ~ Space Optimization

    public class count_subsets_with_sum_K_space_optimized {
        public int findWays(int num[], int k) {
            int n = num.length;

            int prev[] = new int[k + 1];

            prev[0] = 1;

            if (num[0] <= k)
                prev[num[0]] = 1;

            for (int ind = 1; ind < n; ind++) {
                int cur[] = new int[k + 1];
                cur[0] = 1;
                for (int target = 1; target <= k; target++) {

                    int notTaken = prev[target];

                    int taken = 0;
                    if (num[ind] <= target)
                        taken = prev[target - num[ind]];

                    cur[target] = notTaken + taken;
                }
                prev = cur;
            }
            return prev[k];
        }
    }

    // ! Lecture 18

    // ^ Count partitions with given difference

    // ~ Memoization Solution

    public class count_partition_with_given_diff_memoization {
        int mod = 100000007;

        public int countPartitions(int n, int d, int[] arr) {
            int sum = 0;
            for (int i : arr)
                sum += i;
            if (sum - d < 0 || (sum - d) % 2 != 0)
                return 0;
            return findWays(arr, (sum - d) / 2);
        }

        public int findWays(int num[], int tar) {
            int n = num.length;
            int[][] dp = new int[n][tar + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, tar, num, dp);
        }

        public int f(int index, int sum, int[] a, int[][] dp) {
            // * base case written for (Array elements do not include zeros)
            // if (sum == 0)
            // return 1;
            // if (index == 0) {
            // if (a[index] == sum)
            // return 1;
            // else
            // return 0;
            // }
            // * Base case written for Array elements containing zeros
            // this below line creates problem, so we comment it out
            // if (sum == 0)return 1;
            if (index == 0) {
                if (sum == 0 && a[0] == 0)
                    return 2; // 2 possible ways, considering zero and not considering
                if (sum == 0 || sum == a[0])
                    return 1; // only one possible way, not pick
                else
                    return 0;
            }
            if (dp[index][sum] != -1)
                return dp[index][sum];

            int not_pick = f(index - 1, sum, a, dp);
            int pick = 0;
            if (sum >= a[index])
                pick = f(index - 1, sum - a[index], a, dp);
            return dp[index][sum] = (pick + not_pick) % mod;
        }
    }

    // ~ Tabulation

    public class count_partition_with_given_diff_tabulation {
        int mod = (int) (Math.pow(10, 9) + 7);

        public int countPartitions(int n, int d, int[] arr) {
            int totSum = 0;
            for (int i = 0; i < n; i++) {
                totSum += arr[i];
            }

            // Checking for edge cases
            if (totSum - d < 0 || (totSum - d) % 2 == 1)
                return 0;

            return findWays(arr, (totSum - d) / 2);
        }

        public int findWays(int[] num, int tar) {
            int n = num.length;

            int dp[][] = new int[n][tar + 1];

            if (num[0] == 0)
                dp[0][0] = 2; // 2 cases -pick and not pick
            else
                dp[0][0] = 1; // 1 case - not pick

            if (num[0] != 0 && num[0] <= tar)
                dp[0][num[0]] = 1; // 1 case -pick

            for (int ind = 1; ind < n; ind++) {
                for (int target = 0; target <= tar; target++) {

                    int notTaken = dp[ind - 1][target];

                    int taken = 0;
                    if (num[ind] <= target)
                        taken = dp[ind - 1][target - num[ind]];

                    dp[ind][target] = (notTaken + taken) % mod;
                }
            }
            return dp[n - 1][tar];
        }
    }

    // ~ Space Optimized

    public class count_partition_with_given_diff_space_optimized {
        public int mod = (int) (Math.pow(10, 9) + 7);

        public int findWays(int[] num, int tar) {
            int n = num.length;

            int[] prev = new int[tar + 1];

            if (num[0] == 0)
                prev[0] = 2; // 2 cases -pick and not pick
            else
                prev[0] = 1; // 1 case - not pick

            if (num[0] != 0 && num[0] <= tar)
                prev[num[0]] = 1; // 1 case -pick

            for (int ind = 1; ind < n; ind++) {
                int[] curr = new int[tar + 1];
                for (int target = 0; target <= tar; target++) {

                    int notTaken = prev[target];

                    int taken = 0;
                    if (num[ind] <= target)
                        taken = prev[target - num[ind]];

                    curr[target] = (notTaken + taken) % mod;
                }
                prev = curr;
            }
            return prev[tar];
        }

        public int countPartitions(int n, int d, int[] arr) {
            int totSum = 0;
            for (int i = 0; i < n; i++) {
                totSum += arr[i];
            }

            // Checking for edge cases
            if (totSum - d < 0 || (totSum - d) % 2 == 1)
                return 0;

            return findWays(arr, (totSum - d) / 2);
        }
    }

    // ! Lecture 19

    // ^ 0/1 Knapsack

    // A thief wants to rob a store. He is carrying a bag of capacity W. The store
    // has ‘n’ items. Its weight is given by the ‘wt’ array and its value by the
    // ‘val’ array. He can either include an item in its knapsack or exclude it but
    // can’t partially have it as a fraction. We need to find the maximum value of
    // items that the thief can steal.

    // ~ Recursion

    public class knapsack_recursion {
        public int knapsack(int[] weight, int[] value, int n, int maxWeight) {
            return f(n - 1, maxWeight, weight, value);
        }

        public int f(int index, int w, int[] wt, int[] val) {
            if (index == 0) {
                if (w >= wt[0])
                    return val[0];
                else
                    return 0;
            }

            int not_pick = 0 + f(index - 1, w, wt, val);
            int pick = Integer.MIN_VALUE;
            if (w >= wt[index])
                pick = val[index] + f(index - 1, w - wt[index], wt, val);
            return Math.max(pick, not_pick);
        }
    }

    // ~ Memoization

    public class knapsack_memoization {
        public int knapsack(int[] weight, int[] value, int n, int maxWeight) {
            int[][] dp = new int[n][maxWeight + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, maxWeight, weight, value, dp);
        }

        public int f(int index, int w, int[] wt, int[] val, int[][] dp) {
            if (index == 0) {
                if (w >= wt[0])
                    return val[0];
                else
                    return 0;
            }
            if (dp[index][w] != -1)
                return dp[index][w];
            int not_pick = 0 + f(index - 1, w, wt, val, dp);
            int pick = Integer.MIN_VALUE;
            if (w >= wt[index])
                pick = val[index] + f(index - 1, w - wt[index], wt, val, dp);
            return dp[index][w] = Math.max(pick, not_pick);
        }
    }

    // ~ Tabulation

    public class knapsack_tabulation {
        public int knapsack(int[] weight, int[] value, int n, int maxWeight) {
            int[][] dp = new int[n][maxWeight + 1];

            for (int w = weight[0]; w <= maxWeight; w++)
                dp[0][w] = value[0];

            for (int index = 1; index < n; index++) {
                for (int w = 0; w <= maxWeight; w++) {
                    // just copy paste the recurrence below
                    int not_pick = 0 + dp[index - 1][w];
                    int pick = Integer.MIN_VALUE;
                    if (w >= weight[index])
                        pick = value[index] + dp[index - 1][w - weight[index]];
                    dp[index][w] = Math.max(pick, not_pick);
                }
            }
            return dp[n - 1][maxWeight];
        }
    }

    // ~ Space Optimized

    public class knapsack_space_optimized {
        public int knapsack(int[] weight, int[] value, int n, int maxWeight) {
            // int[][] dp = new int[n][maxWeight+1];
            int[] prev = new int[maxWeight + 1];

            for (int w = weight[0]; w <= maxWeight; w++)
                prev[w] = value[0];

            for (int index = 1; index < n; index++) {
                int[] curr = new int[maxWeight + 1];
                for (int w = 0; w <= maxWeight; w++) {
                    // just copy paste the recurrence below
                    int not_pick = 0 + prev[w];
                    int pick = Integer.MIN_VALUE;
                    ;
                    if (w >= weight[index])
                        pick = value[index] + prev[w - weight[index]];
                    curr[w] = Math.max(pick, not_pick);
                }
                prev = curr;
            }
            return prev[maxWeight];
        }
    }

    // ~ Ultimate Space Optimized

    public class Solution {
        public int knapsack(int[] weight, int[] value, int n, int maxWeight) {
            // int[][] dp = new int[n][maxWeight+1];
            int[] prev = new int[maxWeight + 1];

            for (int w = weight[0]; w <= maxWeight; w++)
                prev[w] = value[0];

            for (int index = 1; index < n; index++) {
                // reverse the below loop
                for (int w = maxWeight; w >= 0; w--) {
                    // just copy paste the recurrence below
                    int not_pick = 0 + prev[w];
                    int pick = Integer.MIN_VALUE;
                    if (w >= weight[index])
                        pick = value[index] + prev[w - weight[index]];
                    prev[w] = Math.max(pick, not_pick);
                }
            }
            return prev[maxWeight];
        }
    }

    // ! Lecture 20

    // ^ Minimum Coins

    // We are given a target sum of ‘X’ and ‘N’ distinct numbers denoting the coin
    // denominations. We need to tell the minimum number of coins required to reach
    // the target sum. We can pick a coin denomination for any number of times we
    // want.

    // ~ Recursion

    public class minimum_coins_recursion {
        public int minimumElements(int num[], int x) {
            int n = num.length;
            int ans = f(n - 1, x, num);
            // if ans is not possible
            if (ans >= (int) Math.pow(10, 9))
                return -1;
            else
                return ans;
        }

        public int f(int index, int T, int[] coins) {
            if (index == 0) {
                if (T % coins[0] == 0)
                    return T / coins[0];
                else
                    return (int) Math.pow(10, 9);
            }
            int not_take = 0 + f(index - 1, T, coins);
            int take = (int) Math.pow(10, 9);
            if (T >= coins[index])
                take = 1 + f(index, T - coins[index], coins);

            return Math.min(take, not_take);
        }
    }

    // ~ Memoization

    public class minimum_coins_memoization {
        public int minimumElements(int num[], int x) {
            int n = num.length;
            int[][] dp = new int[n][x + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            int ans = f(n - 1, x, num, dp);
            // if ans is not possible
            if (ans >= (int) Math.pow(10, 9))
                return -1;
            else
                return ans;
        }

        public int f(int index, int T, int[] coins, int[][] dp) {
            if (index == 0) {
                if (T % coins[0] == 0)
                    return T / coins[0];
                else
                    return (int) Math.pow(10, 9);
            }
            if (dp[index][T] != -1)
                return dp[index][T];
            int not_take = 0 + f(index - 1, T, coins, dp);
            int take = (int) Math.pow(10, 9);
            if (T >= coins[index])
                take = 1 + f(index, T - coins[index], coins, dp);

            return dp[index][T] = Math.min(take, not_take);
        }
    }

    // ~ Tabulation

    public class minimum_coins_tabulation {
        public int minimumElements(int coins[], int target) {
            int n = coins.length;
            int[][] dp = new int[n][target + 1];
            // Base Case
            for (int T = 0; T <= target; T++) {
                if (T % coins[0] == 0)
                    dp[0][T] = T / coins[0];
                else
                    dp[0][T] = (int) Math.pow(10, 9);
            }
            for (int index = 1; index < n; index++) {
                for (int T = 0; T <= target; T++) {
                    // copy the recurrence
                    int not_take = 0 + dp[index - 1][T];
                    int take = (int) Math.pow(10, 9);
                    if (T >= coins[index])
                        take = 1 + dp[index][T - coins[index]];
                    dp[index][T] = Math.min(take, not_take);
                }
            }
            int ans = dp[n - 1][target];
            // if ans is not possible
            if (ans >= (int) Math.pow(10, 9))
                return -1;
            else
                return ans;
        }
    }

    // ~ Space Optimized

    public class minimum_coins_space_optimized {
        public int minimumElements(int coins[], int target) {
            int n = coins.length;
            // int[][] dp = new int[n][target+1];
            int[] prev = new int[target + 1];
            int[] curr = new int[target + 1];
            // Base Case
            for (int T = 0; T <= target; T++) {
                if (T % coins[0] == 0)
                    prev[T] = T / coins[0];
                else
                    prev[T] = (int) Math.pow(10, 9);
            }
            for (int index = 1; index < n; index++) {
                for (int T = 0; T <= target; T++) {
                    // copy the recurrence
                    int not_take = 0 + prev[T];
                    int take = (int) Math.pow(10, 9);
                    if (T >= coins[index])
                        take = 1 + curr[T - coins[index]];
                    curr[T] = Math.min(take, not_take);
                }
                prev = curr;
            }
            int ans = prev[target];
            // if ans is not possible
            if (ans >= (int) Math.pow(10, 9))
                return -1;
            else
                return ans;
        }
    }

    // ! Lecture 21

    // ^ Target Sum

    // We are given an array ‘ARR’ of size ‘N’ and a number ‘Target’. Our task is to
    // build an expression from the given array where we can place a ‘+’ or ‘-’ sign
    // in front of an integer. We want to place a sign in front of every integer of
    // the array and get our required target. We need to count the number of ways in
    // which we can achieve our required target.

    public class TargetSum {
        public int targetSum(int n, int target, int[] arr) {
            // & just call the lecture 18 method :)
            return countPartitions(n, target, arr);
        }

        public int findWays(int[] num, int tar) {
            int n = num.length;

            int[] prev = new int[tar + 1];

            if (num[0] == 0)
                prev[0] = 2; // 2 cases -pick and not pick
            else
                prev[0] = 1; // 1 case - not pick

            if (num[0] != 0 && num[0] <= tar)
                prev[num[0]] = 1; // 1 case -pick

            for (int ind = 1; ind < n; ind++) {
                int[] curr = new int[tar + 1];
                for (int target = 0; target <= tar; target++) {

                    int notTaken = prev[target];

                    int taken = 0;
                    if (num[ind] <= target)
                        taken = prev[target - num[ind]];

                    curr[target] = (notTaken + taken);
                }
                prev = curr;
            }
            return prev[tar];
        }

        public int countPartitions(int n, int d, int[] arr) {
            int totSum = 0;
            for (int i = 0; i < n; i++) {
                totSum += arr[i];
            }

            // Checking for edge cases
            if (totSum - d < 0 || (totSum - d) % 2 == 1)
                return 0;

            return findWays(arr, (totSum - d) / 2);
        }
    }

    // ! Lecture 22

    // ^ Coin Change 2

    // We are given an array Arr with N distinct coins and a target. We have an
    // infinite supply of each coin denomination. We need to find the number of ways
    // we sum up the coin values to give us the target.

    // ~ Recursion

    public class coin_change_2_recursion {

        public long countWaysToMakeChange(int denominations[], int value) {
            int n = denominations.length;
            return f(n - 1, value, denominations);
        }

        public long f(int index, int T, int[] a) {
            if (index == 0) {
                return T % a[0] == 0 ? 1 : 0;
            }
            long not_take = 0 + f(index - 1, T, a);
            long take = 0;
            if (T >= a[index])
                take = f(index, T - a[index], a);
            return take + not_take;
        }
    }

    // ~ Memoization

    public class coin_change_2_memoization {

        public long countWaysToMakeChange(int denominations[], int value) {
            int n = denominations.length;
            long[][] dp = new long[n][value + 1];
            for (long[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, value, denominations, dp);
        }

        public long f(int index, int T, int[] a, long[][] dp) {
            if (index == 0) {
                return T % a[0] == 0 ? 1 : 0;
            }
            if (dp[index][T] != -1)
                return dp[index][T];
            long not_take = 0 + f(index - 1, T, a, dp);
            long take = 0;
            if (T >= a[index])
                take = f(index, T - a[index], a, dp);
            return dp[index][T] = take + not_take;
        }
    }

    // ~ Tabulation

    public class coin_change_2_tabulation {

        public long countWaysToMakeChange(int denominations[], int target) {
            int n = denominations.length;
            long[][] dp = new long[n][target + 1];

            for (int T = 0; T <= target; T++) {
                dp[0][T] = T % denominations[0] == 0 ? 1 : 0;
            }

            for (int index = 1; index < n; index++) {
                for (int T = 0; T <= target; T++) {
                    // copy the recurrence
                    long not_take = 0 + dp[index - 1][T];
                    long take = 0;
                    if (T >= denominations[index])
                        take = dp[index][T - denominations[index]];
                    dp[index][T] = take + not_take;
                }
            }
            return dp[n - 1][target];
        }
    }

    // ~ Space Optimized

    public class coin_change_2_space_optimized {

        public long countWaysToMakeChange(int denominations[], int target) {
            int n = denominations.length;
            // long[][] dp = new long[n][target+1];
            long[] prev = new long[target + 1];
            long[] curr = new long[target + 1];

            for (int T = 0; T <= target; T++) {
                prev[T] = T % denominations[0] == 0 ? 1 : 0;
            }

            for (int index = 1; index < n; index++) {
                for (int T = 0; T <= target; T++) {
                    // copy the recurrence
                    long not_take = 0 + prev[T];
                    long take = 0;
                    if (T >= denominations[index])
                        take = curr[T - denominations[index]];
                    curr[T] = take + not_take;
                }
                prev = curr;
            }
            return prev[target];
        }
    }

    // ! Lecture 23

    // ^ Unbound Knapsack

    // A thief wants to rob a store. He is carrying a bag of capacity W. The store
    // has ‘n’ items of infinite supply. Its weight is given by the ‘wt’ array and
    // its value by the ‘val’ array. He can either include an item in its knapsack
    // or exclude it but can’t partially have it as a fraction. We need to find the
    // maximum value of items that the thief can steal. He can take a single item
    // any number of times he wants and put it in his knapsack.

    // ~ Recursion

    public class unbounded_knapsack_recursion {
        public int unboundedKnapsack(int n, int w, int[] profit, int[] weight) {
            return f(n - 1, w, weight, profit);
        }

        public int f(int index, int w, int[] weight, int[] value) {
            if (index == 0)
                return (w / weight[0]) * value[0];
            int not_take = 0 + f(index - 1, w, weight, value);
            int take = (int) Math.pow(-10, 9);
            if (w >= weight[index])
                take = value[index] + f(index, w - weight[index], weight, value);
            return Math.max(take, not_take);
        }
    }

    // ~ Memoization

    public class unbounded_knapsack_memoization {
        public int unboundedKnapsack(int n, int w, int[] profit, int[] weight) {
            int[][] dp = new int[n][w + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, w, weight, profit, dp);
        }

        public int f(int index, int w, int[] weight, int[] value, int[][] dp) {
            if (index == 0)
                return (w / weight[0]) * value[0];
            if (dp[index][w] != -1)
                return dp[index][w];
            int not_take = 0 + f(index - 1, w, weight, value, dp);
            int take = (int) Math.pow(-10, 9);
            if (w >= weight[index])
                take = value[index] + f(index, w - weight[index], weight, value, dp);
            return dp[index][w] = Math.max(take, not_take);
        }
    }

    // ~ Tabulation

    public class unbounded_knapsack_tabulation {
        public int unboundedKnapsack(int n, int wt, int[] value, int[] weight) {
            int[][] dp = new int[n][wt + 1];
            // base case
            for (int w = 0; w <= wt; w++) {
                dp[0][w] = (int) (w / weight[0]) * value[0];
            }
            // writing nested loops according to changing parameters
            for (int index = 1; index < n; index++) {
                for (int w = 0; w <= wt; w++) {
                    // copy the recurrence
                    int not_take = 0 + dp[index - 1][w];
                    int take = (int) Math.pow(-10, 9);
                    if (w >= weight[index])
                        take = value[index] + dp[index][w - weight[index]];
                    dp[index][w] = Math.max(take, not_take);
                }
            }
            return dp[n - 1][wt];
        }
    }

    // ~ Space Optimization

    public class unbounded_knapsack_space_optimization {
        public int unboundedKnapsack(int n, int wt, int[] value, int[] weight) {
            // int[][] dp = new int[n][wt+1];
            int[] prev = new int[wt + 1];
            int[] curr = new int[wt + 1];
            // base case
            for (int w = 0; w <= wt; w++) {
                prev[w] = (int) (w / weight[0]) * value[0];
            }
            // writing nested loops according to changing parameters
            for (int index = 1; index < n; index++) {
                for (int w = 0; w <= wt; w++) {
                    // copy the recurrence
                    int not_take = 0 + prev[w];
                    int take = (int) Math.pow(-10, 9);
                    if (w >= weight[index])
                        take = value[index] + curr[w - weight[index]];
                    curr[w] = Math.max(take, not_take);
                }
                prev = curr;
            }
            return prev[wt];
        }
    }

    // ~ Ultimate Space Optimized

    public class unbounded_knapsack_ultimate_space_optimization {
        public int unboundedKnapsack(int n, int wt, int[] value, int[] weight) {
            // int[][] dp = new int[n][wt+1];
            int[] prev = new int[wt + 1];
            // int[] curr = new int[wt+1];
            // base case
            for (int w = 0; w <= wt; w++) {
                prev[w] = (int) (w / weight[0]) * value[0];
            }
            // writing nested loops according to changing parameters
            for (int index = 1; index < n; index++) {
                for (int w = 0; w <= wt; w++) {
                    // copy the recurrence
                    int not_take = 0 + prev[w];
                    int take = (int) Math.pow(-10, 9);
                    if (w >= weight[index])
                        take = value[index] + prev[w - weight[index]];
                    prev[w] = Math.max(take, not_take);
                }
                // prev=curr;
            }
            return prev[wt];
        }
    }

    // ! Lecture 24

    // ^ Rod Cutting

    // We are given a rod of size ‘N’. It can be cut into pieces. Each length of a
    // piece has a particular price given by the price array. Our task is to find
    // the maximum revenue that can be generated by selling the rod after cutting(
    // if required) into pieces.

    // ~ Recursion

    public class rod_cutting_recursion {
        public int cutRod(int price[], int n) {
            return f(n - 1, n, price);
        }

        public int f(int index, int n, int[] price) {
            if (index == 0)
                return n * price[0];
            int not_take = 0 + f(index - 1, n, price);
            int take = Integer.MIN_VALUE;
            int rod_length = index + 1;
            if (rod_length <= n)
                take = price[index] + f(index, n - rod_length, price);
            return Math.max(take, not_take);
        }
    }

    // ~ Memoization

    public class rod_cutting_memoization {
        public int cutRod(int price[], int n) {
            int[][] dp = new int[n][n + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, n, price, dp);
        }

        public int f(int index, int n, int[] price, int[][] dp) {
            if (index == 0)
                return n * price[0];
            if (dp[index][n] != -1)
                return dp[index][n];
            int not_take = 0 + f(index - 1, n, price, dp);
            int take = Integer.MIN_VALUE;
            int rod_length = index + 1;
            if (rod_length <= n)
                take = price[index] + f(index, n - rod_length, price, dp);
            return dp[index][n] = Math.max(take, not_take);
        }
    }

    // ~ Tabulation

    public class rod_cutting_tabulation {
        public int cutRod(int price[], int size) {
            int[][] dp = new int[size][size + 1];
            // Base case
            for (int i = 0; i <= size; i++)
                dp[0][i] = i * price[0];

            // creating nested loops by checking changing parameters
            for (int index = 1; index < size; index++) {
                for (int n = 0; n <= size; n++) {
                    // copy the recurrence
                    int not_take = 0 + dp[index - 1][n];
                    int take = Integer.MIN_VALUE;
                    int rod_length = index + 1;
                    if (rod_length <= n)
                        take = price[index] + dp[index][n - rod_length];
                    dp[index][n] = Math.max(take, not_take);
                }
            }
            return dp[size - 1][size];
        }
    }

    // ~ Space Optimization

    public class rod_cutting_space_optimized {
        public int cutRod(int price[], int size) {
            // int[][] dp = new int[size][size+1];
            int[] prev = new int[size + 1];
            int[] curr = new int[size + 1];
            // Base case
            for (int i = 0; i <= size; i++)
                prev[i] = i * price[0];

            // creating nested loops by checking changing parameters
            for (int index = 1; index < size; index++) {
                for (int n = 0; n <= size; n++) {
                    // copy the recurrence
                    int not_take = 0 + prev[n];
                    int take = Integer.MIN_VALUE;
                    int rod_length = index + 1;
                    if (rod_length <= n)
                        take = price[index] + curr[n - rod_length];
                    curr[n] = Math.max(take, not_take);
                }
                prev = curr;
            }
            return prev[size];
        }
    }

    // ~ Ultimate Space Optimization

    public class rod_cutting_ultimate_space_optimized {
        public int cutRod(int price[], int size) {
            // int[][] dp = new int[size][size+1];
            int[] prev = new int[size + 1];
            // int[] curr = new int[size+1];
            // Base case
            for (int i = 0; i <= size; i++)
                prev[i] = i * price[0];

            // creating nested loops by checking changing parameters
            for (int index = 1; index < size; index++) {
                for (int n = 0; n <= size; n++) {
                    // copy the recurrence
                    int not_take = 0 + prev[n];
                    int take = Integer.MIN_VALUE;
                    int rod_length = index + 1;
                    if (rod_length <= n)
                        take = price[index] + prev[n - rod_length];
                    prev[n] = Math.max(take, not_take);
                }
                // prev=curr;
            }
            return prev[size];
        }
    }

    // ! Lecture 25

    // ^ Longest Common Subsequence

    // The longest Common Subsequence is defined for two strings. It is the common
    // subsequence that has the greatest length.

    // ~ Recursion

    public class LCS_recursion {

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            return f(m - 1, n - 1, s, t);
        }

        public int f(int i, int j, String s, String t) {
            if (i < 0 || j < 0)
                return 0;
            if (s.charAt(i) == t.charAt(j))
                return 1 + f(i - 1, j - 1, s, t);
            return Math.max(f(i - 1, j, s, t), f(i, j - 1, s, t));
        }
    }

    // ~ Memoization

    public class LCS_memoization {

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[][] dp = new int[m][n];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(m - 1, n - 1, s, t, dp);
        }

        public int f(int i, int j, String s, String t, int[][] dp) {
            if (i < 0 || j < 0)
                return 0;
            if (dp[i][j] != -1)
                return dp[i][j];
            if (s.charAt(i) == t.charAt(j))
                return dp[i][j] = 1 + f(i - 1, j - 1, s, t, dp);
            return dp[i][j] = Math.max(f(i - 1, j, s, t, dp), f(i, j - 1, s, t, dp));
        }
    }

    // ~ Tabulation

    public class LCS_tabulation {

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[][] dp = new int[m + 1][n + 1];

            // base case
            for (int i = 0; i <= m; i++)
                dp[i][0] = 0;
            for (int j = 0; j <= n; j++)
                dp[0][j] = 0;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        dp[i][j] = 1 + dp[i - 1][j - 1];
                    else
                        dp[i][j] = Math.max(dp[i - 1][j], dp[i][j - 1]);
                }
            }

            return dp[m][n];
        }
    }

    // ~ Space Optimization

    public class LCS_Space_optimized {

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[] prev = new int[n + 1];
            int[] curr = new int[n + 1];

            // base case is not needed since array is already initialized to 0

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        curr[j] = 1 + prev[j - 1];
                    else
                        curr[j] = Math.max(prev[j], curr[j - 1]);
                }
                // idk why prev = curr wasn't working
                prev = (int[]) (curr.clone());
            }

            return prev[n];
        }
    }

    // ! Lecture 26

    // ^ Printing LCS

    // ~ Tabulation

    class printing_LCS {

        public void lcs(String s1, String s2) {

            int n = s1.length();
            int m = s2.length();

            int dp[][] = new int[n + 1][m + 1];
            for (int i = 0; i <= n; i++) {
                dp[i][0] = 0;
            }
            for (int i = 0; i <= m; i++) {
                dp[0][i] = 0;
            }

            for (int ind1 = 1; ind1 <= n; ind1++) {
                for (int ind2 = 1; ind2 <= m; ind2++) {
                    if (s1.charAt(ind1 - 1) == s2.charAt(ind2 - 1))
                        dp[ind1][ind2] = 1 + dp[ind1 - 1][ind2 - 1];
                    else
                        dp[ind1][ind2] = 0 + Math.max(dp[ind1 - 1][ind2], dp[ind1][ind2 - 1]);
                }
            }

            int len = dp[n][m];
            int i = n;
            int j = m;

            int index = len - 1;
            String str = "";
            for (int k = 1; k <= len; k++) {
                str += "$"; // dummy string
            }
            StringBuilder ss = new StringBuilder(s1);
            StringBuilder str2 = new StringBuilder(str);
            while (i > 0 && j > 0) {
                if (ss.charAt(i - 1) == s2.charAt(j - 1)) {
                    str2.setCharAt(index, ss.charAt(i - 1));
                    index--;
                    i--;
                    j--;
                } else if (ss.charAt(i - 1) > s2.charAt(j - 1)) {
                    i--;
                } else
                    j--;
            }
            System.out.println(str2);
        }
    }

    // ! Lecture 27

    // ^ Longest Common Substring

    // ~ Tabulation

    public class LC_Substring_tabulation {
        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[][] dp = new int[m + 1][n + 1];

            // base case
            for (int i = 0; i <= m; i++)
                dp[i][0] = 0;
            for (int j = 0; j <= n; j++)
                dp[0][j] = 0;

            int maxi = 0;
            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) { // make i to i-1 and j to j-1
                        dp[i][j] = 1 + dp[i - 1][j - 1];
                        maxi = Math.max(maxi, dp[i][j]);
                    } else
                        dp[i][j] = 0;
                }
            }

            return maxi;
        }
    }

    // ~ Space Optimized

    public class LC_Substring_Space_Optimization {
        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[] prev = new int[n + 1];
            int[] curr = new int[n + 1];

            // base case
            // for(int i = 0; i<=m;i++) dp[i][0] = 0;
            // for(int j = 0; j<=n;j++) prev[j] = 0;

            int maxi = 0;
            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) { // make i to i-1 and j to j-1
                        curr[j] = 1 + prev[j - 1];
                        maxi = Math.max(maxi, curr[j]);
                    } else
                        curr[j] = 0;
                }
                prev = (int[]) (curr.clone());
            }

            return maxi;
        }
    }

    // ! Lecture 28

    // ^ Longest Palindromic Subsequence

    // ~ Space optimization

    public class longest_palindromic_subsequence {
        public int longestPalindromeSubsequence(String s) {
            String ss = new StringBuilder(s).reverse().toString();
            return lcs(s, ss);
        }

        // just the same LCS function of Lecture 25
        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[] prev = new int[n + 1];
            int[] curr = new int[n + 1];

            // base case
            // for(int i = 0; i<=m;i++) dp[i][0] = 0;
            // for(int j = 0; j<=n;j++) prev[j] = 0;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        curr[j] = 1 + prev[j - 1];
                    else
                        curr[j] = Math.max(prev[j], curr[j - 1]);
                }
                prev = (int[]) (curr.clone());
            }

            return prev[n];
        }
    }

    // ! Lecture 29

    // ^ Minimum Insertions to make String Palindrome

    // ~ Space Optimized

    public class min_insertion_to_make_string_palindrome {
        public int minInsertion(String str) {
            int n = str.length();
            int longest_palindromic_subsequence = longestPalindromeSubsequence(str);
            return n - longest_palindromic_subsequence;
        }

        public int longestPalindromeSubsequence(String s) {
            String ss = new StringBuilder(s).reverse().toString();
            return lcs(s, ss);
        }

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[] prev = new int[n + 1];
            int[] curr = new int[n + 1];

            // base case
            // for(int i = 0; i<=m;i++) dp[i][0] = 0;
            // for(int j = 0; j<=n;j++) prev[j] = 0;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        curr[j] = 1 + prev[j - 1];
                    else
                        curr[j] = Math.max(prev[j], curr[j - 1]);
                }
                prev = (int[]) (curr.clone());
            }

            return prev[n];
        }
    }

    // ! Lecture 30

    // ^ Minimum Operations to convert String a to String b

    // ~ Space Optimization

    public class min_op_to_convert_string_a_to_String_b {
        public int canYouMake(String str, String ptr) {
            int str_len = str.length();
            int ptr_len = ptr.length();
            int lcs_length = lcs(str, ptr);

            return str_len + ptr_len - (2 * lcs_length);
        }

        public int lcs(String s, String t) {
            int m = s.length();
            int n = t.length();
            int[] prev = new int[n + 1];
            int[] curr = new int[n + 1];

            // base case
            // for(int i = 0; i<=m;i++) dp[i][0] = 0;
            // for(int j = 0; j<=n;j++) prev[j] = 0;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        curr[j] = 1 + prev[j - 1];
                    else
                        curr[j] = Math.max(prev[j], curr[j - 1]);
                }
                prev = (int[]) (curr.clone());
            }

            return prev[n];
        }
    }

    // ! Lecture 31

    // ^ Shortest Common Super-sequence

    // ~ Tabulation

    public class shortest_Common_Supersequence {
        public String shortestSupersequence(String s, String t) {

            // just the LCS Code
            int m = s.length();
            int n = t.length();
            int[][] dp = new int[m + 1][n + 1];

            // base case
            for (int i = 0; i <= m; i++)
                dp[i][0] = 0;
            for (int j = 0; j <= n; j++)
                dp[0][j] = 0;

            for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) // make i to i-1 and j to j-1
                        dp[i][j] = 1 + dp[i - 1][j - 1];
                    else
                        dp[i][j] = Math.max(dp[i - 1][j], dp[i][j - 1]);
                }
            }

            // new code starts here
            String ans = "";
            int i = m, j = n;

            while (i > 0 && j > 0) {
                if (s.charAt(i - 1) == t.charAt(j - 1)) {
                    ans += s.charAt(i - 1);
                    i--;
                    j--;
                } else if (dp[i - 1][j] > dp[i][j - 1]) {
                    ans += s.charAt(i - 1);
                    i--;
                } else {
                    ans += t.charAt(j - 1);
                    j--;
                }
            }

            // copying the last left over character
            while (i > 0) {
                ans += s.charAt(i - 1);
                i--;
            }
            while (j > 0) {
                ans += t.charAt(j - 1);
                j--;
            }

            String answer = new StringBuilder(ans).reverse().toString();
            return answer;
        }
    }

    // ! Lecture 32

    // ^ Distinct Subsequences

    // We are given two strings S1 and S2, we want to know how many distinct
    // subsequences of S2 are present in S1.

    // ~ Recursion

    public class distinct_subsequences_recursion {
        public int numDistinct(String s, String t) {
            int n = s.length();
            int m = t.length();

            return f(n - 1, m - 1, s, t);
        }

        private int f(int i, int j, String s, String t) {
            if (j < 0)
                return 1; // if t string is exhausted
            if (i < 0)
                return 0; // if s string is exhausted

            if (s.charAt(i) == t.charAt(j)) {
                int leaveOne = f(i - 1, j - 1, s, t); // move both indexes of both strings
                int stay = f(i - 1, j, s, t); // keep the t string index and search for other occurrences
                return leaveOne + stay; // sum up the possibilities
            } else {
                return f(i - 1, j, s, t);
            }
        }
    }

    // ~ Memoization

    public class distinct_subsequences_memoization {
        public int numDistinct(String s, String t) {
            int n = s.length();
            int m = t.length();

            int[][] dp = new int[n][m];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(n - 1, m - 1, s, t, dp);
        }

        private int f(int i, int j, String s, String t, int[][] dp) {
            if (j < 0)
                return 1; // if t string is exhausted
            if (i < 0)
                return 0; // if s string is exhausted

            if (dp[i][j] != -1)
                return dp[i][j];
            if (s.charAt(i) == t.charAt(j)) {
                int leaveOne = f(i - 1, j - 1, s, t, dp); // move both indexes of both strings
                int stay = f(i - 1, j, s, t, dp); // keep the t string index and search for other occurences
                return dp[i][j] = leaveOne + stay; // sum up the possibilities
            } else {
                return dp[i][j] = f(i - 1, j, s, t, dp);
            }
        }
    }

    // ~ Tabulation

    class distinct_subsequences_tabulation {
        public int numDistinct(String s, String t) {
            int n = s.length();
            int m = t.length();

            int[][] dp = new int[n + 1][m + 1];
            // for(int[] row : dp) Arrays.fill(row,-1);

            // base case
            for (int i = 0; i <= n; i++)
                dp[i][0] = 1;
            // 2nd base case line of 'j' is not needed because array is already initialized
            // to 0

            // writing loops according to changing parameters
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) {
                        int leaveOne = dp[i - 1][j - 1]; // move both indexes of both strings
                        int stay = dp[i - 1][j]; // keep the t string index and search for other occurrences
                        dp[i][j] = leaveOne + stay; // sum up the possibilities
                    } else {
                        dp[i][j] = dp[i - 1][j];
                    }
                }
            }
            return dp[n][m];
        }
    }

    // ~ Space Optimization

    class distinct_subsequences_space_optimized {
        public int numDistinct(String s, String t) {
            int n = s.length();
            int m = t.length();

            int[] prev = new int[m + 1];
            int[] curr = new int[m + 1];
            // for(int[] row : dp) Arrays.fill(row,-1);

            // base case
            // for(int i = 0; i<=n; i++) dp[i][0] = 1;
            // 2nd base case line of 'j' is not needed because array is already initialized
            // to 0
            prev[0] = curr[0] = 1;

            // writing loops according to changing parameters
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) {
                        int leaveOne = prev[j - 1]; // move both indexes of both strings
                        int stay = prev[j]; // keep the t string index and search for other occurrences
                        curr[j] = leaveOne + stay; // sum up the possibilities
                    } else {
                        curr[j] = prev[j];
                    }
                }
                prev = curr;
            }
            return prev[m];
        }
    }

    // ~ Ultimate Space Optimized

    class distinct_subsequences_ultimate_space_optimized {
        public int numDistinct(String s, String t) {
            int n = s.length();
            int m = t.length();

            int[] prev = new int[m + 1];
            // for(int[] row : dp) Arrays.fill(row,-1);

            // base case
            // for(int i = 0; i<=n; i++) dp[i][0] = 1;
            // 2nd base case line of 'j' is not needed because array is already initialized
            // to 0
            prev[0] = 1;

            // writing loops according to changing parameters
            for (int i = 1; i <= n; i++) {
                for (int j = m; j >= 1; j--) {
                    // copy the recurrence
                    if (s.charAt(i - 1) == t.charAt(j - 1)) {
                        int leaveOne = prev[j - 1]; // move both indexes of both strings
                        int stay = prev[j]; // keep the t string index and search for other occurrences
                        prev[j] = leaveOne + stay; // sum up the possibilities
                    }
                    // else part is not needed now
                }
            }
            return prev[m];
        }
    }

    // ! Lecture 33

    // ^ Edit Distance

    // We are given two strings ‘S1’ and ‘S2’. We need to convert S1 to S2. The
    // following three operations are allowed:

    // Deletion of a character.
    // Replacement of a character with another one.
    // Insertion of a character.
    // We have to return the minimum number of operations required to convert S1 to
    // S2 as our answer.

    // ~ Recursion

    public class edit_distance_recursion {

        public int editDistance(String str1, String str2) {
            int n = str1.length();
            int m = str2.length();

            return f(n - 1, m - 1, str1, str2);
        }

        public int f(int i, int j, String str1, String str2) {
            if (i < 0)
                return j + 1;
            if (j < 0)
                return i + 1;

            if (str1.charAt(i) == str2.charAt(j)) {
                return 0 + f(i - 1, j - 1, str1, str2);
            } else {
                int insert = 1 + f(i, j - 1, str1, str2);
                int delete = 1 + f(i - 1, j, str1, str2);
                int replace = 1 + f(i - 1, j - 1, str1, str2);

                return Math.min(replace, Math.min(insert, delete));
            }
        }
    }

    // ~ Memoization

    public class edit_distance_memoization {

        public int editDistance(String str1, String str2) {
            int n = str1.length();
            int m = str2.length();
            int[][] dp = new int[n][m];
            for (int[] row : dp)
                Arrays.fill(row, -1);

            return f(n - 1, m - 1, str1, str2, dp);
        }

        public int f(int i, int j, String str1, String str2, int[][] dp) {
            if (i < 0)
                return j + 1;
            if (j < 0)
                return i + 1;

            if (dp[i][j] != -1)
                return dp[i][j];
            if (str1.charAt(i) == str2.charAt(j)) {
                return dp[i][j] = 0 + f(i - 1, j - 1, str1, str2, dp);
            } else {
                int insert = 1 + f(i, j - 1, str1, str2, dp);
                int delete = 1 + f(i - 1, j, str1, str2, dp);
                int replace = 1 + f(i - 1, j - 1, str1, str2, dp);

                return dp[i][j] = Math.min(replace, Math.min(insert, delete));
            }
        }
    }

    // ~ Tabulation

    public class edit_distance_tabulation {

        public int editDistance(String str1, String str2) {
            int n = str1.length();
            int m = str2.length();
            int[][] dp = new int[n + 1][m + 1];
            // for(int[] row : dp) Arrays.fill(row,-1);

            for (int i = 0; i <= n; i++)
                dp[i][0] = i;
            for (int j = 0; j <= m; j++)
                dp[0][j] = j;

            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if (str1.charAt(i - 1) == str2.charAt(j - 1)) {
                        dp[i][j] = dp[i - 1][j - 1];
                    } else {
                        int insert = 1 + dp[i][j - 1];
                        int delete = 1 + dp[i - 1][j];
                        int replace = 1 + dp[i - 1][j - 1];

                        dp[i][j] = Math.min(replace, Math.min(insert, delete));
                    }
                }
            }

            return dp[n][m];
        }
    }

    // ~ Space Optimization

    public class edit_distance_space_optimized {

        public int editDistance(String str1, String str2) {
            int n = str1.length();
            int m = str2.length();
            int[] prev = new int[m + 1];
            int[] curr = new int[m + 1];
            // for(int[] row : dp) Arrays.fill(row,-1);

            // for(int i = 0; i<=n; i++) dp[i][0] = i;
            for (int j = 0; j <= m; j++)
                prev[j] = j;

            for (int i = 1; i <= n; i++) {
                curr[0] = i;
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if (str1.charAt(i - 1) == str2.charAt(j - 1)) {
                        curr[j] = prev[j - 1];
                    } else {
                        int insert = 1 + curr[j - 1];
                        int delete = 1 + prev[j];
                        int replace = 1 + prev[j - 1];

                        curr[j] = Math.min(replace, Math.min(insert, delete));
                    }
                }
                // idk why we need to do this
                prev = (int[]) (curr.clone());
            }

            return prev[m];
        }
    }

    // ! Lecture 34

    // ^ Wildcard Matching

    // We are given two strings ‘S1’ and ‘S2’. String S1 can have the following two
    // special characters:

    // ‘?’ can be matched to a single character of S2.
    // ‘*’ can be matched to any sequence of characters of S2. (sequence can be of
    // length zero or more).
    // We need to check whether strings S1 and S2 match or not.

    // ~ Recursion

    public class wildcard_matching_recursion {
        public boolean wildcardMatching(String pattern, String text) {
            int n = pattern.length();
            int m = text.length();
            return f(n - 1, m - 1, pattern, text);
        }

        public boolean f(int i, int j, String pattern, String text) {
            // base case
            if (i < 0 && j < 0)
                return true;
            if (i < 0 && j >= 0)
                return false;
            if (j < 0 && i >= 0) {
                for (int ii = 0; ii <= i; ii++) {
                    if (pattern.charAt(ii) != '*')
                        return false;
                }
                return true;
            }
            // trying all possible ways
            if ((pattern.charAt(i) == text.charAt(j)) || pattern.charAt(i) == '?')
                return f(i - 1, j - 1, pattern, text);
            else if (pattern.charAt(i) == '*')
                return f(i - 1, j, pattern, text) | f(i, j - 1, pattern, text);
            else
                return false;
        }
    }

    // ~ Memoization

    public class wildcard_matching_memoization {
        public boolean wildcardMatching(String pattern, String text) {
            int n = pattern.length();
            int m = text.length();
            // Boolean is different from boolean
            Boolean[][] dp = new Boolean[n][m];
            return f(n - 1, m - 1, pattern, text, dp);
        }

        public boolean f(int i, int j, String pattern, String text, Boolean[][] dp) {
            if (i < 0 && j < 0)
                return true;
            if (i < 0 && j >= 0)
                return false;
            if (j < 0 && i >= 0) {
                for (int ii = 0; ii <= i; ii++) {
                    if (pattern.charAt(ii) != '*')
                        return false;
                }
                return true;
            }

            if (dp[i][j] != null)
                return dp[i][j];
            if ((pattern.charAt(i) == text.charAt(j)) || pattern.charAt(i) == '?')
                return dp[i][j] = f(i - 1, j - 1, pattern, text, dp);
            else if (pattern.charAt(i) == '*')
                return dp[i][j] = f(i - 1, j, pattern, text, dp) | f(i, j - 1, pattern, text, dp);
            else
                return dp[i][j] = false;
        }
    }

    // ~ Memoization in 1-Based Indexing

    public class wildcard_matching_memoization_1_based {
        public boolean wildcardMatching(String pattern, String text) {
            int n = pattern.length();
            int m = text.length();
            Boolean[][] dp = new Boolean[n + 1][m + 1];
            return f(n, m, pattern, text, dp);
        }

        public boolean f(int i, int j, String pattern, String text, Boolean[][] dp) {
            if (i == 0 && j == 0)
                return true;
            if (i == 0 && j > 0)
                return false;
            if (j == 0 && i > 0) {
                for (int ii = 1; ii <= i; ii++) {
                    if (pattern.charAt(ii - 1) != '*')
                        return false;
                }
                return true;
            }

            if (dp[i][j] != null)
                return dp[i][j];
            if ((pattern.charAt(i - 1) == text.charAt(j - 1)) || pattern.charAt(i - 1) == '?')
                return dp[i][j] = f(i - 1, j - 1, pattern, text, dp);
            else if (pattern.charAt(i - 1) == '*')
                return dp[i][j] = f(i - 1, j, pattern, text, dp) | f(i, j - 1, pattern, text, dp);
            else
                return dp[i][j] = false;
        }
    }

    // ~ Tabulation

    public class wildcard_matching_tabulation {
        public boolean wildcardMatching(String pattern, String text) {
            int n = pattern.length();
            int m = text.length();
            Boolean[][] dp = new Boolean[n + 1][m + 1];

            // base case
            dp[0][0] = true;
            for (int j = 1; j <= m; j++)
                dp[0][j] = false;
            for (int i = 0; i <= n; i++) {
                boolean flag = true;
                for (int ii = 1; ii <= i; ii++) {
                    if (pattern.charAt(ii - 1) != '*') {
                        flag = false;
                        break;
                    }
                }
                dp[i][0] = flag;
            }
            // make nested loops
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if ((pattern.charAt(i - 1) == text.charAt(j - 1)) || pattern.charAt(i - 1) == '?')
                        dp[i][j] = dp[i - 1][j - 1];
                    else if (pattern.charAt(i - 1) == '*')
                        dp[i][j] = dp[i - 1][j] | dp[i][j - 1];
                    else
                        dp[i][j] = false;
                }
            }

            return dp[n][m];
        }
    }

    // ~ Space Optimization

    public class wildcard_matching_space_optimized {
        public boolean wildcardMatching(String pattern, String text) {
            int n = pattern.length();
            int m = text.length();
            Boolean[] prev = new Boolean[m + 1];
            Boolean[] curr = new Boolean[m + 1];
            // prev is always the 0th row, and curr is the next row
            // base case
            prev[0] = true;
            for (int j = 1; j <= m; j++)
                prev[j] = false;

            // make nested loops
            for (int i = 1; i <= n; i++) {
                curr[0] = isAllStars(pattern, i);
                for (int j = 1; j <= m; j++) {
                    // copy the recurrence
                    if ((pattern.charAt(i - 1) == text.charAt(j - 1)) || pattern.charAt(i - 1) == '?')
                        curr[j] = prev[j - 1];
                    else if (pattern.charAt(i - 1) == '*')
                        curr[j] = prev[j] || curr[j - 1];
                    else
                        curr[j] = false;
                }
                // idk why prev = curr won't work
                prev = (Boolean[]) (curr.clone());
            }

            return prev[m];
        }

        public boolean isAllStars(String S1, int i) {

            // S1 is taken in 1-based indexing
            for (int j = 1; j <= i; j++) {
                if (S1.charAt(j - 1) != '*')
                    return false;
            }
            return true;
        }
    }

    // ! Lecture 35

    // ^ Best Time to buy and Sell Stocks 1

    // ~ Normal Iterative Solution

    public class bet_time_to_buy_and_sell_stocks {
        public int maximumProfit(ArrayList<Integer> prices) {
            int mini = Integer.MAX_VALUE;
            int profit = 0;
            for (int i : prices) {
                if (i <= mini) {
                    mini = i;
                } else {
                    profit = Math.max(profit, i - mini);
                }
            }
            return profit;
        }
    }

    // ! Lecture 36

    // ^ Best Time to Buy and Sell Stocks 2

    // We are given an array Arr[] of length n. It represents the price of a stock
    // on ‘n’ days. The following guidelines need to be followed:

    // We can buy and sell the stock any number of times.
    // In order to sell the stock, we need to first buy it on the same or any
    // previous day.
    // We can’t buy a stock again after buying it once. In other words, we first buy
    // a stock and then sell it. After selling we can buy and sell again. But we
    // can’t sell before buying and can’t buy before selling any previously bought
    // stock.

    // ~ Recursion

    public class bet_time_to_buy_and_sell_stocks_2_recursion {
        public long getMaximumProfit(int n, long[] prices) {
            return f(0, 1, prices, n);
        }

        public long f(int index, int buy, long[] prices, int n) {
            // base case
            if (index == n)
                return 0;

            // if we can buy (we don't hold a stock)
            if (buy == 1) {
                long will_buy = -prices[index] + f(index + 1, 0, prices, n);
                long not_buy = 0 + f(index + 1, 1, prices, n);
                return Math.max(will_buy, not_buy);
            }
            // if we cannot buy (we are holding a stock)
            else {
                long will_sell = +prices[index] + f(index + 1, 1, prices, n);
                long not_sell = 0 + f(index + 1, 0, prices, n);
                return Math.max(will_sell, not_sell);
            }
        }
    }

    // ~ Memoization

    public class bet_time_to_buy_and_sell_stocks_2_memoization {
        public long getMaximumProfit(int n, long[] prices) {
            long[][] dp = new long[n][2];
            for (long[] row : dp)
                Arrays.fill(row, -1);
            return f(0, 1, prices, n, dp);
        }

        public long f(int index, int buy, long[] prices, int n, long[][] dp) {
            // base case
            if (index == n)
                return 0;

            if (dp[index][buy] != -1)
                return dp[index][buy];

            long profit = 0;
            // if we can buy (we don't hold a stock)
            if (buy == 1) {
                long will_buy = -prices[index] + f(index + 1, 0, prices, n, dp);
                long not_buy = 0 + f(index + 1, 1, prices, n, dp);
                profit = Math.max(will_buy, not_buy);
            }
            // if we cannot buy (we are holding a stock)
            else {
                long will_sell = +prices[index] + f(index + 1, 1, prices, n, dp);
                long not_sell = 0 + f(index + 1, 0, prices, n, dp);
                profit = Math.max(will_sell, not_sell);
            }
            return dp[index][buy] = profit;
        }
    }

    // ~ Tabulation

    public class bet_time_to_buy_and_sell_stocks_2_tabulation {
        public long getMaximumProfit(int n, long[] prices) {
            long[][] dp = new long[n + 1][2];
            // Base Case
            dp[n][0] = dp[n][1] = 0;

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int buy = 0; buy <= 1; buy++) {
                    // copy the recurrence
                    long profit = 0;
                    // if we can buy (we don't hold a stock)
                    if (buy == 1) {
                        long will_buy = -prices[index] + dp[index + 1][0];
                        long not_buy = 0 + dp[index + 1][1];
                        profit = Math.max(will_buy, not_buy);
                    }
                    // if we cannot buy (we are holding a stock)
                    else {
                        long will_sell = +prices[index] + dp[index + 1][1];
                        long not_sell = 0 + dp[index + 1][0];
                        profit = Math.max(will_sell, not_sell);
                    }
                    dp[index][buy] = profit;
                }
            }
            return dp[0][1];
        }
    }

    // ~ Space Optimization

    public class bet_time_to_buy_and_sell_stocks_2_space_optimization {
        public long getMaximumProfit(int n, long[] prices) {
            long[] ahead = new long[2];
            long[] curr = new long[2];
            // Base Case
            ahead[0] = ahead[1] = 0;

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int buy = 0; buy <= 1; buy++) {
                    // copy the recurrence
                    long profit = 0;
                    // if we can buy (we don't hold a stock)
                    if (buy == 1) {
                        long will_buy = -prices[index] + ahead[0];
                        long not_buy = 0 + ahead[1];
                        profit = Math.max(will_buy, not_buy);
                    }
                    // if we cannot buy (we are holding a stock)
                    else {
                        long will_sell = +prices[index] + ahead[1];
                        long not_sell = 0 + ahead[0];
                        profit = Math.max(will_sell, not_sell);
                    }
                    curr[buy] = profit;
                }
                ahead = curr;
            }
            return ahead[1];
        }
    }

    // ! Lecture 37

    // ^ Best Time to Buy and Sell Stocks 3

    // We are given an array Arr[] of length n. It represents the price of a stock
    // on ‘n’ days. The following guidelines need to be followed:

    // We can buy and sell the stock any number of times.
    // In order to sell the stock, we need to first buy it on the same or any
    // previous day.
    // We can’t buy a stock again after buying it once. In other words, we first buy
    // a stock and then sell it. After selling we can buy and sell again. But we
    // can’t sell before buying and can’t buy before selling any previously bought
    // stock.
    // We can do at most 2 transactions.

    // ~ Recursion

    public class bet_time_to_buy_and_sell_stocks_3_recursion {
        public int maxProfit(ArrayList<Integer> prices, int n) {
            return f(0, 1, 2, prices, n);
        }

        public int f(int index, int can_buy, int cap, List<Integer> prices, int n) {
            // base case
            if (index == n)
                return 0;
            if (cap == 0)
                return 0;

            int profit = 0;
            // try all ways
            if (can_buy == 1) {
                int will_buy = -prices.get(index) + f(index + 1, 0, cap, prices, n);
                int not_buy = 0 + f(index + 1, 1, cap, prices, n);
                profit = Math.max(will_buy, not_buy);
            } else {
                int will_sell = prices.get(index) + f(index + 1, 1, cap - 1, prices, n);
                int not_sell = 0 + f(index + 1, 0, cap, prices, n);
                profit = Math.max(will_sell, not_sell);
            }
            return profit;
        }
    }

    // ~ Memoization

    public class bet_time_to_buy_and_sell_stocks_3_memoization {
        public int maxProfit(ArrayList<Integer> prices, int n) {
            int[][][] dp = new int[n][2][3];
            for (int[][] row2 : dp) {
                for (int[] row1 : row2)
                    Arrays.fill(row1, -1);
            }
            return f(0, 1, 2, prices, n, dp);
        }

        public int f(int index, int can_buy, int cap, List<Integer> prices,
                int n, int[][][] dp) {
            // base case
            if (index == n)
                return 0;
            if (cap == 0)
                return 0;

            if (dp[index][can_buy][cap] != -1)
                return dp[index][can_buy][cap];

            int profit = 0;
            // try all ways
            if (can_buy == 1) {
                int will_buy = -prices.get(index) + f(index + 1, 0, cap, prices, n, dp);
                int not_buy = 0 + f(index + 1, 1, cap, prices, n, dp);
                profit = Math.max(will_buy, not_buy);
            } else {
                int will_sell = prices.get(index) + f(index + 1, 1, cap - 1, prices, n, dp);
                int not_sell = 0 + f(index + 1, 0, cap, prices, n, dp);
                profit = Math.max(will_sell, not_sell);
            }
            return dp[index][can_buy][cap] = profit;
        }
    }

    // ~ Tabulation

    public class bet_time_to_buy_and_sell_stocks_3_tabulation {
        public int maxProfit(ArrayList<Integer> prices, int n) {
            int[][][] dp = new int[n + 1][2][3];

            // writing base cases won't make sense because we have already
            // initialized the array to 0;

            // writing nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int can_buy = 0; can_buy <= 1; can_buy++) {
                    // we ignore cap = 0 because it is already covered in base case
                    for (int cap = 1; cap <= 2; cap++) {
                        // copy the recurrence
                        int profit = 0;
                        // try all ways
                        if (can_buy == 1) {
                            int will_buy = -prices.get(index) + dp[index + 1][0][cap];
                            int not_buy = 0 + dp[index + 1][1][cap];
                            profit = Math.max(will_buy, not_buy);
                        } else {
                            int will_sell = prices.get(index) + dp[index + 1][1][cap - 1];
                            int not_sell = 0 + dp[index + 1][0][cap];
                            profit = Math.max(will_sell, not_sell);
                        }
                        dp[index][can_buy][cap] = profit;
                    }
                }
            }
            return dp[0][1][2];
        }
    }

    // ~ Space Optimization

    public class bet_time_to_buy_and_sell_stocks_3_space_optimization {
        public int maxProfit(ArrayList<Integer> prices, int n) {
            int[][] ahead = new int[2][3];
            int[][] curr = new int[2][3];
            // replace dp[index+1] with ahead

            // writing base cases won't make sense because we have already
            // initialized the array to 0;

            // writing nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int can_buy = 0; can_buy <= 1; can_buy++) {
                    // we ignore cap = 0 because it is already covered in base case
                    for (int cap = 1; cap <= 2; cap++) {
                        // copy the recurrence
                        int profit = 0;
                        // try all ways
                        if (can_buy == 1) {
                            int will_buy = -prices.get(index) + ahead[0][cap];
                            int not_buy = 0 + ahead[1][cap];
                            profit = Math.max(will_buy, not_buy);
                        } else {
                            int will_sell = prices.get(index) + ahead[1][cap - 1];
                            int not_sell = 0 + ahead[0][cap];
                            profit = Math.max(will_sell, not_sell);
                        }
                        curr[can_buy][cap] = profit;
                    }
                }
                ahead = curr;
            }
            return ahead[1][2];
        }
    }

    // ! Lecture 38

    // ^ Best Time to Buy and Sell Stocks 4

    // Same Previous Question, just instead of Limitation of 2 Transactions, we are
    // Limited to K Transactions

    // ~ Space Optimized

    public class bet_time_to_buy_and_sell_stocks_4_space_optimization {
        public int maximumProfit(int[] prices, int n, int k) {
            int[][] ahead = new int[2][k + 1];
            int[][] curr = new int[2][k + 1];
            // replace dp[index+1] with ahead

            // writing base cases won't make sense because we have already
            // initialized the array to 0;

            // writing nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int can_buy = 0; can_buy <= 1; can_buy++) {
                    // we ignore cap = 0 because it is already covered in base case
                    for (int cap = 1; cap <= k; cap++) {
                        // copy the recurrence
                        int profit = 0;
                        // try all ways
                        if (can_buy == 1) {
                            int will_buy = -prices[index] + ahead[0][cap];
                            int not_buy = 0 + ahead[1][cap];
                            profit = Math.max(will_buy, not_buy);
                        } else {
                            int will_sell = prices[index] + ahead[1][cap - 1];
                            int not_sell = 0 + ahead[0][cap];
                            profit = Math.max(will_sell, not_sell);
                        }
                        curr[can_buy][cap] = profit;
                    }
                }
                ahead = curr;
            }
            return ahead[1][k];
        }
    }

    // ! Lecture 39

    // ^ Best Time to Buy and Sell Stocks with Cool-down

    // just the same Buy and Sell Stocks 2 question but we cannot buy just after we
    // sold, we need to wait for one day

    // ~ Tabulation

    public class bet_time_to_buy_and_sell_stocks_with_CoolDown {

        public int stockProfit(int[] prices) {
            int n = prices.length;
            // create n+2 size array instead of n+1
            // because n+2 index is visited
            int[][] dp = new int[n + 2][2];
            // Base Case
            dp[n][0] = dp[n][1] = 0;

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int buy = 0; buy <= 1; buy++) {
                    // copy the recurrence
                    int profit = 0;
                    // if we can buy (we don't hold a stock)
                    if (buy == 1) {
                        int will_buy = -prices[index] + dp[index + 1][0];
                        int not_buy = 0 + dp[index + 1][1];
                        profit = Math.max(will_buy, not_buy);
                    }
                    // if we cannot buy (we are holding a stock)
                    else {
                        // int will_sell = +prices[index] + dp[index + 1][1];
                        // one small change only from Buy and sell stocks 2 i.e. (index+2)
                        int will_sell = +prices[index] + dp[index + 2][1];
                        int not_sell = 0 + dp[index + 1][0];
                        profit = Math.max(will_sell, not_sell);
                    }
                    dp[index][buy] = profit;
                }
            }
            return dp[0][1];
        }
    }
    // Cannot be Space optimized because we are dealing with index+1 and index+2
    // so hypothetically we would need 3 arrays, which is not feasible

    // ~ Another Approach by removing the 'buy' wala for loop

    public class bet_time_to_buy_and_sell_stocks_with_CoolDown_One_less_loop {

        public int stockProfit(int[] prices) {
            int n = prices.length;
            // create n+2 size array instead of n+1
            // because n+2 index is visited
            int[][] dp = new int[n + 2][2];
            // Base Case
            dp[n][0] = dp[n][1] = 0;

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {

                // copy the recurrence
                // we can omit the 'buy' wala loop and remove the if and else conditions

                // when buy = 1, dp[index][1] will be executed

                // if we can buy (we don't hold a stock)

                int will_buy = -prices[index] + dp[index + 1][0];
                int not_buy = 0 + dp[index + 1][1];
                dp[index][1] = Math.max(will_buy, not_buy);

                // if we cannot buy (we are holding a stock)

                // when buy == 0, dp[index][0] will be executed
                int will_sell = +prices[index] + dp[index + 2][1];
                int not_sell = 0 + dp[index + 1][0];
                dp[index][0] = Math.max(will_sell, not_sell);

            }

            return dp[0][1];
        }
    }

    // Now this can be Space Optimized

    // ~ Space Optimization

    public class bet_time_to_buy_and_sell_stocks_with_CoolDown_space_optimized {

        public int stockProfit(int[] prices) {
            int n = prices.length;
            // create n+2 size array instead of n+1
            // because n+2 index is visited
            int[] front2 = new int[2];
            int[] front1 = new int[2];
            int[] curr = new int[2];
            // Base Case
            // dp[n][0] = dp[n][1] = 0;

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {

                // we can omit the 'buy' wala loop and remove the if and else conditions

                // when buy == 1, dp[index][1] will be executed

                // if we can buy (we don't hold a stock)

                int will_buy = -prices[index] + front1[0];
                int not_buy = 0 + front1[1];
                curr[1] = Math.max(will_buy, not_buy);

                // if we cannot buy (we are holding a stock)

                // when buy == 0, dp[index][0] will be executed
                int will_sell = +prices[index] + front2[1];
                int not_sell = 0 + front1[0];
                curr[0] = Math.max(will_sell, not_sell);

                front2 = (int[]) (front1.clone());
                front1 = (int[]) (curr.clone());

            }

            return curr[1];
        }
    }

    // ! Lecture 40

    // ^ Best Time to Buy and Sell Stocks with Transaction Fees

    // ~ Space Optimized

    public class bet_time_to_buy_and_sell_stocks_with_transaction_fees {
        public int maximumProfit(int n, int fees, int[] prices) {
            int[] ahead = new int[2];
            int[] curr = new int[2];
            // Base Case not needed because array is already initialized to 0

            // write nested loops
            for (int index = n - 1; index >= 0; index--) {

                // we can omit the 'buy' wala loop and remove the if and else conditions

                // when buy == 1, dp[index][1] will be executed

                // if we can buy (we don't hold a stock)

                int will_buy = -prices[index] + ahead[0];
                int not_buy = 0 + ahead[1];
                curr[1] = Math.max(will_buy, not_buy);

                // if we cannot buy (we are holding a stock)

                // when buy == 0, dp[index][0] will be executed
                int will_sell = +prices[index] - fees + ahead[1];
                int not_sell = 0 + ahead[0];
                curr[0] = Math.max(will_sell, not_sell);

                ahead = (int[]) (curr.clone());

            }

            return curr[1];
        }
    }

    // ! Lecture 41

    // ^ Longest Increasing Subsequence

    // ~ Recursion

    // will give TLE
    public class LIS_recursion {

        public int longestIncreasingSubsequence(int arr[]) {
            return f(0, -1, arr);
        }

        public int f(int index, int prev, int[] arr) {
            int n = arr.length;
            int len;
            if (index == n)
                return 0;

            // not take
            len = 0 + f(index + 1, prev, arr);
            // take
            if (prev == -1 || arr[index] > arr[prev]) {
                len = Math.max(len, 1 + f(index + 1, index, arr));
            }

            return len;
        }
    }

    // ~ Memoization

    // will give Runtime Error because of the constraints, we have to make
    // 10^5 x 10^5 size grid which is not possible
    public class LIS_memoization {
        public int longestIncreasingSubsequence(int arr[]) {
            int n = arr.length;
            int[][] dp = new int[n][n + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(0, -1, arr, dp);
        }

        public int f(int index, int prev, int[] arr, int[][] dp) {
            int n = arr.length;
            int len;
            if (index == n)
                return 0;

            if (dp[index][prev + 1] != -1)
                return dp[index][prev + 1];
            // not take
            len = 0 + f(index + 1, prev, arr, dp);
            // take
            if (prev == -1 || arr[index] > arr[prev]) {
                len = Math.max(len, 1 + f(index + 1, index, arr, dp));
            }

            return dp[index][prev + 1] = len;
        }
    }

    // ! Lecture 42

    // ^ Printing LIS and Tabulation of LIS

    // ~ Tabulation

    // will give Runtime Error because of the constraints, we have to make
    // 10^5 x 10^5 size grid which is not possible
    public class LIS_Tabulation {
        public int longestIncreasingSubsequence(int arr[]) {
            int n = arr.length, len;
            int[][] dp = new int[n + 1][n + 1];

            // make nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int prev = index - 1; prev >= -1; prev--) {
                    // copy the recurrence

                    // not take
                    len = 0 + dp[index + 1][prev + 1];
                    // take
                    if (prev == -1 || arr[index] > arr[prev]) {
                        len = Math.max(len, 1 + dp[index + 1][index + 1]);
                    }

                    dp[index][prev + 1] = len;
                }
            }
            return dp[0][-1 + 1];
        }
    }

    // ~ Space Optimization

    // will still give TLE
    public class LIS_space_optimized {
        public int longestIncreasingSubsequence(int arr[]) {
            int n = arr.length, len;
            int[] next = new int[n + 1];
            int[] curr = new int[n + 1];

            // make nested loops
            for (int index = n - 1; index >= 0; index--) {
                for (int prev = index - 1; prev >= -1; prev--) {
                    // copy the recurrence

                    // not take
                    len = 0 + next[prev + 1];
                    // take
                    if (prev == -1 || arr[index] > arr[prev]) {
                        len = Math.max(len, 1 + next[index + 1]);
                    }

                    curr[prev + 1] = len;
                }
                next = curr;
            }
            return next[-1 + 1];
        }
    }

    // ~ Other Solution

    // will still give TLE because of nested loops
    public class LIS_other {
        public int longestIncreasingSubsequence(int arr[]) {
            int n = arr.length, len;
            int[] dp = new int[n];
            Arrays.fill(dp, 1);
            int maxi = 1;
            // make nested loops
            for (int i = 0; i < n; i++) {
                for (int prev = 0; prev < i; prev++) {
                    if (arr[i] > arr[prev])
                        dp[i] = Math.max(dp[i], 1 + dp[prev]);
                }
                maxi = Math.max(maxi, dp[i]);
            }
            return maxi;
        }
    }

    // ! Lecture 44

    // ^ Longest Divisible Subset

    // Pre Requisites - Printing an LIS

    public class LDS {
        public ArrayList<Integer> divisibleSet(int arr[]) {
            int n = arr.length;

            // sort the array

            Arrays.sort(arr);

            int[] dp = new int[n];
            int[] hash = new int[n];
            Arrays.fill(dp, 1);
            Arrays.fill(hash, 1);

            for (int i = 0; i <= n - 1; i++) {

                hash[i] = i; // initializing with current index
                for (int prev_index = 0; prev_index <= i - 1; prev_index++) {

                    if (arr[i] % arr[prev_index] == 0 && 1 + dp[prev_index] > dp[i]) {
                        dp[i] = 1 + dp[prev_index];
                        hash[i] = prev_index;
                    }
                }
            }

            int ans = -1;
            int lastIndex = -1;

            for (int i = 0; i <= n - 1; i++) {
                if (dp[i] > ans) {
                    ans = dp[i];
                    lastIndex = i;
                }
            }

            ArrayList<Integer> temp = new ArrayList<>();
            temp.add(arr[lastIndex]);

            while (hash[lastIndex] != lastIndex) { // till not reach the initialization value
                lastIndex = hash[lastIndex];
                temp.add(arr[lastIndex]);
            }

            // reverse the List
            Collections.reverse(temp);

            return temp;
        }

    }

    // ! Lecture 45

    // ^ Longest String Chain

    // just modifying the LIS function a bit

    class Longest_String_chain {
        public int longestStrChain(int n, String[] arr) {
            // sort according to length
            Arrays.sort(arr, (a, b) -> Integer.compare(a.length(), b.length()));
            // vector<int> dp(n,1);
            int[] dp = new int[n];
            Arrays.fill(dp, 1);

            int maxi = 1;

            for (int i = 0; i <= n - 1; i++) {

                for (int prev_index = 0; prev_index <= i - 1; prev_index++) {

                    if (compare(arr[i], arr[prev_index])) {
                        dp[i] = Math.max(dp[i], 1 + dp[prev_index]);
                    }
                }

                maxi = Math.max(maxi, dp[i]);
            }
            return maxi;
        }

        public boolean compare(String s1, String s2) {
            if (s1.length() != s2.length() + 1)
                return false;

            int first = 0;
            int second = 0;

            while (first < s1.length()) {
                if (second < s2.length() && s1.charAt(first) == s2.charAt(second)) {
                    first++;
                    second++;
                } else
                    first++;
            }
            if (first == s1.length() && second == s2.length())
                return true;
            else
                return false;
        }
    }

    // ! Lecture 46

    // ^ Longest Bi-Tonic Subsequence

    public class Bitonic {
        public int longestBitonicSequence(int arr[], int n) {
            int[] dp1 = new int[n];
            Arrays.fill(dp1, 1);
            int[] dp2 = new int[n];
            Arrays.fill(dp2, 1);

            // LIS dp1 array
            for (int i = 0; i < n; i++) {
                for (int prev = 0; prev < i; prev++) {
                    if (arr[i] > arr[prev])
                        dp1[i] = Math.max(dp1[i], 1 + dp1[prev]);
                }
            }

            // LDS dp2 array
            for (int i = 0; i < n; i++) {
                for (int prev = 0; prev < i; prev++) {
                    // just reverse the sign
                    if (arr[i] < arr[prev])
                        dp1[i] = Math.max(dp1[i], 1 + dp1[prev]);
                }
            }

            // finding bitonic
            int maxi = 0;
            for (int i = 0; i < n; i++) {
                maxi = Math.max(maxi, dp1[i] + dp2[i] - 1);
            }
            return maxi;
        }
    }

    // ! lecture 47

    // ^ Count the number of LIS

    public class count_LIS {
        public int findNumberOfLIS(int n, int arr[]) {
            int[] dp = new int[n];
            Arrays.fill(dp, 1);
            int[] cnt = new int[n];
            Arrays.fill(cnt, 1);
            int maxi = 1;
            // make nested loops
            for (int i = 0; i < n; i++) {
                for (int prev = 0; prev < i; prev++) {
                    if (arr[i] > arr[prev] && 1 + dp[prev] > dp[i]) {
                        dp[i] = 1 + dp[prev];
                        cnt[i] = cnt[prev];
                    } else if (arr[i] > arr[prev] && 1 + dp[prev] == dp[i])
                        cnt[i] += cnt[prev];
                }
                maxi = Math.max(maxi, dp[i]);
            }
            int nos = 0;
            for (int i = 0; i < n; i++) {
                if (dp[i] == maxi)
                    nos += cnt[i];
            }
            return nos;
        }
    }

    // ! Lecture 48

    // ^ Matrix Chain Multiplication

    // Given a chain of matrices, you have to find the minimum cost to multiply
    // these matrices

    // ~ Recursion

    public class MCM_Recursion {

        public int matrixMultiplication(int[] arr, int N) {
            return f(1, N - 1, arr);
        }

        public int f(int i, int j, int[] arr) {
            if (i == j)
                return 0;

            int mini = Integer.MAX_VALUE;
            for (int k = i; k < j; k++) {
                // doing partition and trying out all partitions to find the best
                int steps = (arr[i - 1] * arr[k] * arr[j]) + f(i, k, arr) + f(k + 1, j, arr);
                mini = Math.min(mini, steps);
            }
            return mini;
        }
    }

    // ~ Memoization

    public class MCM_memoization {
        public int matrixMultiplication(int[] arr, int N) {
            int[][] dp = new int[N][N];
            for (int[] row : dp)
                Arrays.fill(row, -1);
            return f(1, N - 1, arr, dp);
        }

        public int f(int i, int j, int[] arr, int[][] dp) {
            if (i == j)
                return 0;

            if (dp[i][j] != -1)
                return dp[i][j];
            int mini = Integer.MAX_VALUE;
            for (int k = i; k < j; k++) {
                int steps = (arr[i - 1] * arr[k] * arr[j]) + f(i, k, arr, dp) + f(k + 1, j, arr, dp);
                mini = Math.min(mini, steps);
            }
            return dp[i][j] = mini;
        }
    }

    // ! Lecture 49

    // ^ Matrix Chain Multiplication Tabulation

    // ~ Tabulation

    public class MCM_Tabulation {
        public int matrixMultiplication(int[] arr, int N) {
            int[][] dp = new int[N][N];

            // base case (can also skip this because dp is already initialized to 0)
            for (int i = 0; i < N; i++)
                dp[i][i] = 0;

            // write nested loops according to changing parameters
            for (int i = N - 1; i >= 1; i--) {
                // we start j from i+1 instead of 1 because i and j are at 2 ends
                // and when they meet, that's the base case
                for (int j = i + 1; j < N; j++) {
                    // copy the recurrence
                    int mini = Integer.MAX_VALUE;
                    for (int k = i; k < j; k++) {
                        int steps = (arr[i - 1] * arr[k] * arr[j]) + dp[i][k] + dp[k + 1][j];
                        mini = Math.min(mini, steps);
                    }
                    dp[i][j] = mini;

                }
            }
            return dp[1][N - 1];
        }
    }

    // ! Lecture 50

    // ^ Minimum Cost to Cut the Stick

    // Given a wooden stick of length n units. The stick is labelled from 0 to n.
    // For example, a stick of length 6 is labelled as follows:

    // Given an integer array cuts where cuts[i] denotes a position you should
    // perform a cut at.

    // You should perform the cuts in order, you can change the order of the cuts as
    // you wish.

    // The cost of one cut is the length of the stick to be cut, the total cost is
    // the sum of costs of all cuts. When you cut a stick, it will be split into two
    // smaller sticks (i.e. the sum of their lengths is the length of the stick
    // before the cut). Please refer to the first example for a better explanation.

    // Return the minimum total cost of the cuts.

    // ~ Recursion

    public class min_cost_to_cut_stick_recursion {
        public int cost(int n, int c, int cuts[]) {
            Arrays.sort(cuts);
            int[] cutz = new int[c + 2];
            cutz[0] = 0; // adding 0 to the starting of the array
            // adding rest of thee elements
            int j = 1;
            for (int i : cuts) {
                cutz[j++] = i;
            }
            cutz[cutz.length - 1] = n; // adding the length at the end

            return f(1, c, cutz);
        }

        public int f(int i, int j, int[] cutz) {
            if (i > j)
                return 0;

            int mini = Integer.MAX_VALUE;
            for (int index = i; index <= j; index++) {
                int cost = cutz[j + 1] - cutz[i - 1] + f(i, index - 1, cutz) + f(index + 1, j, cutz);
                mini = Math.min(mini, cost);
            }
            return mini;
        }
    }

    // ~ Memoization

    public class min_cost_to_cut_stick__memoization {
        public int cost(int n, int c, int cuts[]) {
            Arrays.sort(cuts);
            int[] cutz = new int[c + 2];
            cutz[0] = 0; // adding 0 to the starting of the array
            // adding rest of thee elements
            int j = 1;
            for (int i : cuts) {
                cutz[j++] = i;
            }
            cutz[cutz.length - 1] = n; // adding the length at the end
            int[][] dp = new int[c + 1][c + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);

            return f(1, c, cutz, dp);
        }

        public int f(int i, int j, int[] cutz, int[][] dp) {
            if (i > j)
                return 0;

            if (dp[i][j] != -1)
                return dp[i][j];
            int mini = Integer.MAX_VALUE;
            for (int index = i; index <= j; index++) {
                int cost = cutz[j + 1] - cutz[i - 1] + f(i, index - 1, cutz, dp) + f(index + 1, j, cutz, dp);
                mini = Math.min(mini, cost);
            }
            return dp[i][j] = mini;
        }
    }

    // ~ Tabulation

    public class min_cost_to_cut_stick_tabulation {
        public int cost(int n, int c, int cuts[]) {
            Arrays.sort(cuts);
            int[] cutz = new int[c + 2];
            cutz[0] = 0; // adding 0 to the starting of the array
            // adding rest of thee elements
            int z = 1;
            for (int i : cuts) {
                cutz[z++] = i;
            }
            cutz[cutz.length - 1] = n; // adding the length at the end
            int[][] dp = new int[c + 2][c + 2];
            for (int i = c; i >= 1; i--) {
                for (int j = 1; j <= c; j++) {
                    // copy the recurrence
                    if (i > j)
                        continue;// base case
                    int mini = Integer.MAX_VALUE;
                    for (int index = i; index <= j; index++) {
                        int cost = cutz[j + 1] - cutz[i - 1] + dp[i][index - 1] + dp[index + 1][j];
                        mini = Math.min(mini, cost);
                    }
                    dp[i][j] = mini;
                }
            }
            return dp[1][c];
        }
    }

    // ! Lecture 51

    // ^ Burst Balloons

    // ou are given n balloons, indexed from 0 to n - 1. Each balloon is painted
    // with a number on it represented by an array nums. You are asked to burst all
    // the balloons.

    // If you burst the ith balloon, you will get nums[i - 1] * nums[i] * nums[i +
    // 1] coins. If i - 1 or i + 1 goes out of bounds of the array, then treat it as
    // if there is a balloon with a 1 painted on it.

    // Return the maximum coins you can collect by bursting the balloons wisely.

    // ~ Recursion

    public class burst_balloons_recursion {
        public int maxCoins(int arr[]) {
            int n = arr.length;

            int[] a = new int[n + 2];
            a[0] = 1;
            int z = 1;
            for (int i : arr) {
                a[z++] = i;
            }
            a[a.length - 1] = 1;

            return f(1, n, a);
        }

        public int f(int i, int j, int[] a) {
            if (i > j)
                return 0;

            int maxi = Integer.MIN_VALUE;
            for (int index = i; index <= j; index++) {
                int coins = (a[i - 1] * a[index] * a[j + 1]) + f(i, index - 1, a) + f(index + 1, j, a);
                maxi = Math.max(maxi, coins);
            }
            return maxi;
        }
    }

    // ~ Memoization

    public class burst_balloons_memoization {
        public int maxCoins(int arr[]) {
            int n = arr.length;
            // adding 1 to both ends
            int[] a = new int[n + 2];
            a[0] = 1;
            int z = 1;
            for (int i : arr) {
                a[z++] = i;
            }
            a[a.length - 1] = 1;

            int[][] dp = new int[n + 1][n + 1];
            for (int[] row : dp)
                Arrays.fill(row, -1);

            return f(1, n, a, dp);
        }

        public int f(int i, int j, int[] a, int[][] dp) {
            if (i > j)
                return 0;

            if (dp[i][j] != -1)
                return dp[i][j];
            int maxi = Integer.MIN_VALUE;
            for (int index = i; index <= j; index++) {
                int coins = (a[i - 1] * a[index] * a[j + 1]) + f(i, index - 1, a, dp) + f(index + 1, j, a, dp);
                maxi = Math.max(maxi, coins);
            }
            return dp[i][j] = maxi;
        }
    }

    // ~ Tabulation

    public class burst_balloons_tabulation {
        public int maxCoins(int arr[]) {
            int n = arr.length;
            // adding 1 to both ends
            int[] a = new int[n + 2];
            a[0] = 1;
            int z = 1;
            for (int i : arr) {
                a[z++] = i;
            }
            a[a.length - 1] = 1;

            int[][] dp = new int[n + 2][n + 2]; // make it n+2 for tabulation
            // writing nested loops
            for (int i = n; i >= 1; i--) {
                for (int j = 1; j <= n; j++) {
                    // copy the recurrence
                    if (i > j)
                        continue; // base case
                    int maxi = Integer.MIN_VALUE;
                    for (int index = i; index <= j; index++) {
                        int coins = (a[i - 1] * a[index] * a[j + 1]) + dp[i][index - 1] + dp[index + 1][j];
                        maxi = Math.max(maxi, coins);
                    }
                    dp[i][j] = maxi;
                }
            }

            return dp[1][n];
        }
    }

    // ! Lecture 52

    // ^ Evaluate Boolean Expression to True

    // Find number of ways in which the given boolean expression can be evaluated a
    // true

    // ~ Recursion

    public class boolean_expression_to_true_recursion {
        int mod = 1000000007;

        public int evaluateExp(String exp) {
            int n = exp.length();
            return f(0, n - 1, 1, exp);
        }

        public int f(int i, int j, int isTrue, String exp) {
            if (i > j)
                return 0;
            if (i == j) {
                // when we want true
                if (isTrue == 1)
                    return (exp.charAt(i) == 'T') ? 1 : 0;
                // when we want false
                else
                    return (exp.charAt(i) == 'F') ? 1 : 0;
            }

            long ways = 0;
            for (int index = i + 1; index <= j - 1; index += 2) {
                int LT = f(i, index - 1, 1, exp);
                int LF = f(i, index - 1, 0, exp);
                int RT = f(index + 1, j, 1, exp);
                int RF = f(index + 1, j, 0, exp);

                if (exp.charAt(index) == '&') {
                    if (isTrue == 1) {
                        ways += (LT * RT);
                    } else {
                        ways += (LF * RF) + (LT * RF) + (LF * RT);
                    }
                } else if (exp.charAt(index) == '|') {
                    if (isTrue == 1) {
                        ways += (LT * RT) + (LT * RF) + (LF * RT);
                    } else {
                        ways += (LF * RF);
                    }
                } else {
                    if (isTrue == 1) {
                        ways += (LT * RF) + (LF * RT);
                    } else {
                        ways += (LT * RT) + (LF * RF);
                    }
                }
            }
            return (int) ways % mod;
        }
    }

    // ~ Memoization

    public class boolean_expression_to_true_memoization {
        int mod = 1000000007;

        public int evaluateExp(String exp) {
            int n = exp.length();
            int[][][] dp = new int[n][n][2];
            for (int[][] row1 : dp) {
                for (int[] row2 : row1)
                    Arrays.fill(row2, -1);
            }
            return f(0, n - 1, 1, exp, dp);
        }

        public int f(int i, int j, int isTrue, String exp, int[][][] dp) {
            if (i > j)
                return 0;
            if (i == j) {
                // when we want true
                if (isTrue == 1)
                    return (exp.charAt(i) == 'T') ? 1 : 0;
                // when we want false
                else
                    return (exp.charAt(i) == 'F') ? 1 : 0;
            }

            if (dp[i][j][isTrue] != -1)
                return dp[i][j][isTrue];

            long ways = 0;
            for (int index = i + 1; index <= j - 1; index += 2) {
                int LT = f(i, index - 1, 1, exp, dp);
                int LF = f(i, index - 1, 0, exp, dp);
                int RT = f(index + 1, j, 1, exp, dp);
                int RF = f(index + 1, j, 0, exp, dp);

                if (exp.charAt(index) == '&') {
                    if (isTrue == 1) {
                        ways += (LT * RT);
                    } else {
                        ways += (LF * RF) + (LT * RF) + (LF * RT);
                    }
                } else if (exp.charAt(index) == '|') {
                    if (isTrue == 1) {
                        ways += (LT * RT) + (LT * RF) + (LF * RT);
                    } else {
                        ways += (LF * RF);
                    }
                } else {
                    if (isTrue == 1) {
                        ways += (LT * RF) + (LF * RT);
                    } else {
                        ways += (LT * RT) + (LF * RF);
                    }
                }
            }
            return dp[i][j][isTrue] = (int) ways % mod;
        }
    }

    // ! Lecture 53

    // ^ Palindrome Partitioning 2

    // Given a string s, partition s such that every substring of the partition is a
    // palindrome.

    // Return the minimum cuts needed for a palindrome partitioning of s.

    // ~ Recursion

    public class palindrome_Partitioning_2_recursion {

        public int palindromePartitioning(String str) {
            return f(0, str) - 1;
        }

        public int f(int i, String str) {
            int n = str.length();
            if (i == n)
                return 0;

            String temp = "";
            int mini = Integer.MAX_VALUE;
            for (int j = i; j < n; j++) {
                temp += str.charAt(j);
                if (isPal(temp)) {
                    int cost = 1 + f(j + 1, str);
                    mini = Math.min(mini, cost);
                }
            }
            return mini;
        }

        public boolean isPal(String s) {
            String t = new StringBuilder(s).reverse().toString();
            return s.equals(t);
        }
    }

    // ~ Memoization

    public class palindrome_Partitioning_2_memoization {

        public int palindromePartitioning(String str) {
            int[] dp = new int[str.length()];
            Arrays.fill(dp, -1);
            return f(0, str, dp) - 1;
        }

        public int f(int i, String str, int[] dp) {
            int n = str.length();
            if (i == n)
                return 0;

            if (dp[i] != -1)
                return dp[i];
            String temp = "";
            int mini = Integer.MAX_VALUE;
            for (int j = i; j < n; j++) {
                temp += str.charAt(j);
                if (isPal(temp)) {
                    int cost = 1 + f(j + 1, str, dp);
                    mini = Math.min(mini, cost);
                }
            }
            return dp[i] = mini;
        }

        public boolean isPal(String s) {
            String t = new StringBuilder(s).reverse().toString();
            return s.equals(t);
        }
    }

    // ~ Tabulation

    public class palindrome_Partitioning_2_tabulation {

        public int palindromePartitioning(String str) {
            int n = str.length();
            int[] dp = new int[n + 1];

            for (int i = n - 1; i >= 0; i--) {
                // copy the recurrence
                String temp = "";
                int mini = Integer.MAX_VALUE;
                for (int j = i; j < n; j++) {
                    temp += str.charAt(j);
                    if (isPal(temp)) {
                        int cost = 1 + dp[j + 1];
                        mini = Math.min(mini, cost);
                    }
                }
                dp[i] = mini;
            }
            return dp[0] - 1;
        }

        public boolean isPal(String s) {
            String t = new StringBuilder(s).reverse().toString();
            return s.equals(t);
        }
    }

    // ! Lecture 54

    // ^ Partition Array for Maximum Sum

    // Given an integer array arr, partition the array into (contiguous) subarrays
    // of length at most k. After partitioning, each subarray has their values
    // changed to become the maximum value of that subarray.

    // Return the largest sum of the given array after partitioning.

    // ~ Recursion

    public class partition_array_for_max_sum_recursion {
        public int maximumSubarray(int num[], int k) {
            return f(0, k, num);
        }

        public int f(int index, int k, int[] num) {
            int n = num.length;
            if (index == n)
                return 0;

            int len = 0, maxi = Integer.MIN_VALUE;
            int ans = Integer.MIN_VALUE;
            for (int j = index; j < Math.min(n, index + k); j++) {
                len++;
                maxi = Math.max(maxi, num[j]);
                int sum = (len * maxi) + f(j + 1, k, num);
                ans = Math.max(ans, sum);
            }
            return ans;
        }
    }

    // ~ Memoization

    public class partition_array_for_max_sum_memoization {
        public int maximumSubarray(int num[], int k) {
            int[] dp = new int[num.length];
            Arrays.fill(dp, -1);
            return f(0, k, num, dp);
        }

        public int f(int index, int k, int[] num, int[] dp) {
            int n = num.length;
            if (index == n)
                return 0;

            if (dp[index] != -1)
                return dp[index];
            int len = 0, maxi = Integer.MIN_VALUE;
            int ans = Integer.MIN_VALUE;
            for (int j = index; j < Math.min(n, index + k); j++) {
                len++;
                maxi = Math.max(maxi, num[j]);
                int sum = (len * maxi) + f(j + 1, k, num, dp);
                ans = Math.max(ans, sum);
            }
            return dp[index] = ans;
        }
    }

    // ~ Tabulation

    public class partition_array_for_max_sum_tabulation {
        public int maximumSubarray(int num[], int k) {
            int n = num.length;
            int[] dp = new int[n + 1];

            for (int index = n - 1; index >= 0; index--) {
                // copy the recurrence
                int len = 0, maxi = Integer.MIN_VALUE;
                int ans = Integer.MIN_VALUE;
                for (int j = index; j < Math.min(n, index + k); j++) {
                    len++;
                    maxi = Math.max(maxi, num[j]);
                    int sum = (len * maxi) + dp[j + 1];
                    ans = Math.max(ans, sum);
                }
                dp[index] = ans;
            }
            return dp[0];
        }
    }
}