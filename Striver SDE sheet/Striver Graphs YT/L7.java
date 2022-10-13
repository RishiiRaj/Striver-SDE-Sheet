// ! Lecture 7

// ^ Number of Provinces

// Given an undirected graph with V vertices. We say two vertices u and v belong
// to a single province if there is a path from u to v or v to u. Your task is
// to find the number of provinces.

// Note: A province is a group of directly or indirectly connected cities and no
// other cities outside of the group.

//{ Driver Code Starts
import java.io.*;
import java.util.*;

class L7 {
    public static void main(String args[]) throws IOException {
        BufferedReader read = new BufferedReader(new InputStreamReader(System.in));
        int t = Integer.parseInt(read.readLine());
        while (t-- > 0) {
            int V = Integer.parseInt(read.readLine());

            ArrayList<ArrayList<Integer>> adj = new ArrayList<>();

            for (int i = 0; i < V; i++) {
                String S[] = read.readLine().split(" ");
                ArrayList<Integer> temp = new ArrayList<>();
                for (int j = 0; j < V; j++)
                    temp.add(Integer.parseInt(S[j]));
                adj.add(temp);
            }

            Solution ob = new Solution();
            System.out.println(ob.numProvinces(adj, V));
        }
    }
}
// } Driver Code Ends

// User function Template for Java

class Solution {
    static int numProvinces(ArrayList<ArrayList<Integer>> adj, int V) {
        ArrayList<ArrayList<Integer>> adjLst = new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < V; i++) {
            adjLst.add(new ArrayList<>());
        }

        // changing Adj Matrix to Adj List
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (adj.get(i).get(j) == 1 && i != j) {
                    adjLst.get(i).add(j);
                    adjLst.get(j).add(i);
                }
            }
        }
        // logic starts here

        boolean vis[] = new boolean[V];
        int cnt = 0;
        for (int i = 0; i < V; i++) {
            // we will enter the if condition for each province
            if (!vis[i]) {
                cnt++;
                dfs(i, adjLst, vis);
            }
        }
        return cnt;
    }

    // dfs traversal
    static void dfs(int node, ArrayList<ArrayList<Integer>> adjLst, boolean[] vis) {
        vis[node] = true;
        for (Integer i : adjLst.get(node)) {
            if (!vis[i])
                dfs(i, adjLst, vis);
        }
    }
}