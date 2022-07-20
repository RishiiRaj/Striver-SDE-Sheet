// ! Striver SDE sheet, all Greedy Algorithm solutions

import java.util.*;

class striver_greedy {

    // ? Day 8 Greedy Algorithm

    // ^ N meetings in one room

    // There is one meeting room in a firm. There are N meetings in the form of
    // (start[i], end[i]) where start[i] is start time of meeting i and end[i] is
    // finish time of meeting i.
    // What is the maximum number of meetings that can be accommodated in the
    // meeting room when only one meeting can be held in the meeting room at a
    // particular time?

    // Note: Start time of one chosen meeting can't be equal to the end time of the
    // other chosen meeting.

    // Example 1:

    // Input:
    // N = 6
    // start[] = {1,3,0,5,8,5}
    // end[] = {2,4,6,7,9,9}
    // Output:
    // 4
    // Explanation:
    // Maximum four meetings can be held with
    // given start and end timings.
    // The meetings are - (1, 2),(3, 4), (5,7) and (8,9)

    class Pair {
        int start;
        int end;

        Pair(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }

    class max_meetings {
        public int maxMeetings(int start[], int end[], int n) {
            // add 1st meeting to our answer
            int count = 1;

            Pair arr[] = new Pair[n];
            // convert start[] and end[] arr to pairs
            for (int i = 0; i < n; i++) {
                arr[i] = new Pair(start[i], end[i]);
            }

            // sort pairs according to the end time;
            Arrays.sort(arr, new Comparator<Pair>() {
                public int compare(Pair p1, Pair p2) {
                    return p1.end - p2.end;
                }
            });

            int recentEnded = arr[0].end;

            // if start of current meeting is greater than recentEnded meeting then only add
            // it to our answer
            for (int i = 1; i < n; i++) {
                if (arr[i].start > recentEnded) {
                    count++;
                    recentEnded = arr[i].end;
                }
            }

            return count;
        }
    }

    // ^ Minimum number of platforms required for a railway

    // Given arrival and departure times of all trains that reach a railway station.
    // Find the minimum number of platforms required for the railway station so that
    // no train is kept waiting.
    // Consider that all the trains arrive on the same day and leave on the same
    // day. Arrival and departure time can never be the same for a train but we can
    // have arrival time of one train equal to departure time of the other. At any
    // given instance of time, same platform can not be used for both departure of a
    // train and arrival of another train. In such cases, we need different
    // platforms.

    // Example 1:

    // Input: n = 6
    // arr[] = {0900, 0940, 0950, 1100, 1500, 1800}
    // dep[] = {0910, 1200, 1120, 1130, 1900, 2000}
    // Output: 3
    // Explanation:
    // Minimum 3 platforms are required to
    // safely arrive and depart all trains.

    class max_platforms {
        // Function to find the minimum number of platforms required at the
        // railway station such that no train waits.
        public int findPlatform(int arr[], int dep[], int n) {
            Arrays.sort(arr);
            Arrays.sort(dep);

            // initially platform needed will be 1 for 0 indexed train
            int plat_needed = 1, max = 1;

            // since we already considered the first train, so pointer for
            // arrival will start from 1, and for departure, will start from 0
            int i = 1, j = 0;

            while (i < n && j < n) {
                // we use <= because if dep. time of a train and arrival time of another train
                // is the same, then they will need 2 platforms
                // if arrival time of next train is less than or equal to the departure of the
                // next train
                // a platform will be needed, and i pointer will be incremented
                if (arr[i] <= dep[j]) {
                    plat_needed++;
                    i++;
                }
                // else platform will be vacant, so counter will decrease
                // and departure pointer will be incremented
                else {
                    plat_needed--;
                    j++;
                }
                max = Math.max(max, plat_needed);
            }
            return max;
        }

    }

    // ^ Job sequencing Problem

    // with job if and only if the job is completed by its deadline.

    // Find the number of jobs done and the maximum profit.

    // Note: Jobs will be given in the form (Jobid, Deadline, Profit) associated
    // with that Job.

    // Example 1:

    // Input:
    // N = 4
    // Jobs = {(1,4,20),(2,1,10),(3,1,40),(4,1,30)}
    // Output:
    // 2 60
    // Explanation:
    // Job1 and Job3 can be done with
    // maximum profit of 60 (20+40).

    class Job {
        int id, profit, deadline;

        Job(int id, int profit, int deadline) {
            this.id = id;
            this.profit = profit;
            this.deadline = deadline;
        }
    }

    class job_sequencing {
        // Function to find the maximum profit and the number of jobs done.
        int[] JobScheduling(Job arr[], int n) {
            Arrays.sort(arr, (a, b) -> (b.profit - a.profit));
            int jobs = 0, max_profit = 0;

            // finding maximum deadline
            int maxi = 0;
            for (int i = 0; i < n; i++) {
                if (arr[i].deadline > maxi)
                    maxi = arr[i].deadline;
            }

            int result[] = new int[maxi + 1];

            Arrays.fill(result, -1);

            for (int i = 0; i < n; i++) {
                // this loop is in reverse order because we tend to complete the job
                // with last deadline, at the last, in order to do more jobs
                // in the remaining time, to maximize profit
                for (int j = arr[i].deadline; j > 0; j--) {
                    // if free slot found
                    if (result[j] == -1) {
                        result[j] = i;
                        jobs++;
                        max_profit += arr[i].profit;
                        break;
                    }
                    // if not found, j decreases, and we do that job one day before
                }
            }
            int[] ans = new int[2];
            ans[0] = jobs;
            ans[1] = max_profit;

            return ans;
        }
    }

    // ^ Fractional Knapsack Problem

    // Given weights and values of N items, we need to put these items in a knapsack
    // of capacity W to get the maximum total value in the knapsack.
    // Note: Unlike 0/1 knapsack, you are allowed to break the item.

    // Example 1:

    // Input:
    // N = 3, W = 50
    // values[] = {60,100,120}
    // weight[] = {10,20,30}
    // Output:
    // 220.00
    // Explanation:Total maximum value of item
    // we can have is 220.00 from the given
    // capacity of sack.

    class Item {
        int value, weight;

        Item(int x, int y) {
            this.value = x;
            this.weight = y;
        }
    }

    // comparator class to sort in descending order of value/weight
    class itemCompare implements Comparator<Item> {
        @Override
        public int compare(Item a, Item b) {
            double r1 = (double) (a.value) / (double) (a.weight);
            double r2 = (double) (b.value) / (double) (b.weight);

            if (r1 < r2)
                return 1;
            else if (r1 > r2)
                return -1;
            else
                return 0;
        }
    }

    class fractional_Knapsack {
        // Function to get the maximum total value in the knapsack.
        double fractionalKnapsack(int W, Item arr[], int n) {
            // sorting according to descending order of value/weight ratio
            Arrays.sort(arr, new itemCompare());

            int currWeight = 0;
            double finalValue = 0.0;

            for (int i = 0; i < n; i++) {
                // if current knapsack weight + weight of current item
                // is less than or equal to the capacity of knapsack
                if (currWeight + arr[i].weight <= W) {
                    currWeight += arr[i].weight; // weight is added
                    finalValue += arr[i].value; // it's value is added
                }
                // if not
                else {
                    int remain = W - currWeight; /// we find out the remaining amount of weight
                    finalValue += (double) arr[i].value / arr[i].weight * (double) remain; // and that fraction of value
                                                                                           // is added
                    break;
                }
            }
            return finalValue;
        }
    }

    // ^ Greedy algorithm to find minimum number of coins (CodingNinjas)

    // Dora, the explorer, visits India and decides to try the famous Indian food.
    // However, the restaurant accepts only Indian currency i.e. [1, 2, 5, 10, 20,
    // 50, 100, 500, 1000] valued coins.
    // So, Dora goes to a bank that has an infinite supply of each of the
    // denominations to make a change for a given ‘Amount’ of money. As a cashier at
    // the bank, your task is to provide Dora the minimum number of coins that add
    // up to the given ‘Amount’.

    // For Example
    // For Amount = 70, the minimum number of coins required is 2 i.e an Rs. 50 coin
    // and a Rs. 20 coin.

    public class find_Minimum_Coins {
        public int findMinimumCoins(int amount) {
            ArrayList<Integer> ans = new ArrayList<>();
            int[] coins = { 1, 2, 5, 10, 20, 50, 100, 500, 1000 };
            int len = coins.length;

            for (int i = len - 1; i >= 0; i--) {
                // since coin value can't exceed the amount value
                while (amount >= coins[i]) {
                    amount -= coins[i]; // amount is decreased by that coin value
                    ans.add(coins[i]); // and that coin is added to the list
                }
            }
            return ans.size();
        }
    }

    // ^ Activity Selection (it is the same as N meeting in one room)
}