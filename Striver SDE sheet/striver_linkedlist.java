//! Striver SDE sheet all Linked List Solutions

import java.util.*;

public class striver_linkedlist {

    // ~ Linked List Definition

    public class ListNode {
        int val;
        ListNode next;

        ListNode() {
        }

        ListNode(int val) {
            this.val = val;
        }

        ListNode(int val, ListNode next) {
            this.val = val;
            this.next = next;
        }
    }

    // ? Day 5 Linked List

    // ^ Reverse a Linked List
    // Given the head of a singly linked list, reverse the list, and return the
    // reversed list.
    // Input: head = [1,2,3,4,5]
    // Output: [5,4,3,2,1]

    // 3 pointer approach is used here

    class reverse_ll {
        public ListNode reverseList(ListNode head) {
            if (head == null || head.next == null) {
                return head;
            }
            ListNode prev = head, curr = head.next;
            while (curr != null) {
                ListNode next = curr.next;
                curr.next = prev;
                prev = curr;
                curr = next;
            }
            head.next = null;
            head = prev;
            return head;
        }
    }

    // ^ Find middle of Linked List
    // Given the head of a singly linked list, return the middle node of the linked
    // list.
    // If there are two middle nodes, return the second middle node.

    class middle_of_LL {
        public ListNode middleNode(ListNode head) {
            // using tortoise method
            if (head == null || head.next == null) {
                return head;
            }
            ListNode fast = head, slow = head;
            while (fast != null && fast.next != null) {
                slow = slow.next;
                fast = fast.next.next;

            }
            return slow;
        }
    }

    // ^ Merge 2 Sorted Lists
    // You are given the heads of two sorted linked lists list1 and list2.
    // Merge the two lists in a one sorted list. The list should be made by splicing
    // together the nodes of the first two lists.
    // Return the head of the merged linked list.

    // Input: list1 = [1,2,4], list2 = [1,3,4]
    // Output: [1,1,2,3,4,4]

    class merge_2_sorted_lists {
        public ListNode mergeTwoLists(ListNode l1, ListNode l2) {
            if (l1 == null && l2 == null) {
                return null;
            }
            if (l1 == null) {
                return l2;
            }
            if (l2 == null) {
                return l1;
            }
            ListNode result; // dummy node to store the resultant.
            if (l1.val < l2.val) { // agar l1 ka first element l2 se chhota hai toh usko daalenge result me
                result = l1;
                result.next = mergeTwoLists(l1.next, l2); // fir same function call kiya jayega
                // lekin l1 ka second element aayega ab
            } else { // warna l2 ko daala jayega
                result = l2;
                result.next = mergeTwoLists(l1, l2.next);
            }
            return result;
        }
    }

    // ^ Remove N'th Node from back of the Linked List
    // Given the head of a linked list, remove the nth node from the end of the list
    // and return its head.

    // brute force approach
    class remove_nth_node_from_back_brute {
        public ListNode removeNthFromEnd(ListNode head, int n) {
            if (head.next == null && n == 1 || head == null) {
                return null;
            }
            ListNode currNode = head;

            // count number of nodes
            int count = 1;
            while (currNode.next != null) {
                count += 1;
                currNode = currNode.next;
            }

            int pos = count - n; // the position to be deleted.
            // then we want to remove pos + 1 position
            ListNode curr = head;
            // if it is the first node
            if (pos == 0) {
                head = curr.next;
            } else {
                // else traverse till the required node.
                while (pos - 1 > 0) {
                    curr = curr.next;
                    pos--;
                }
                curr.next = curr.next.next;
            }

            return head;
        }
    }

    // optimal solution

    class remove_nth_node_from_back_optimal {
        public ListNode removeNthFromEnd(ListNode head, int n) {
            if (head == null || head.next == null && n == 1)
                return null;

            ListNode dummy = new ListNode(); // creating a dummy node
            dummy.next = head; // that points to head
            ListNode fast = dummy, slow = dummy; // 2 more pointers pointing to that dummy node
            for (int i = 0; i < n; i++) { // move fast pointer n spaces
                fast = fast.next;
            }
            // base case - if node to be deleted is the last node
            if (fast.next == null)
                return head.next;

            // if that's not the case, move both slo and fast pointer one step till fast
            // reaches the last node
            while (fast.next != null) { // till fast reaches last node
                fast = fast.next;
                slow = slow.next;
            }
            // then unlink the next of slow pointer node.
            slow.next = slow.next.next;
            return head;
        }
    }

    // ^ Add 2 numbers as Linked List
    // You are given two non-empty linked lists representing two non-negative
    // integers. The digits are stored in reverse order, and each of their nodes
    // contains a single digit. Add the two numbers and return the sum as a linked
    // list.
    // You may assume the two numbers do not contain any leading zero, except the
    // number 0 itself.

    class add_2_numbers_as_linked_list {

        public ListNode addTwoNumbers(ListNode l1, ListNode l2) {
            ListNode dummy = new ListNode(); // a dummy node is created
            ListNode temp = dummy; // a temp node is created to iterate
            int carry = 0; // by default the carry is set to zero.

            // while loop will run until l1 is not null, or l2 is not null or carry has some
            // value in it
            while (l1 != null || l2 != null || carry == 1) { // we traverse until l1 is null, l2 is null and carry is
                                                             // zero
                int sum = 0; // by default sum is zero

                // if l1 is not null, its value is added to sum
                if (l1 != null) {
                    sum += l1.val;
                    l1 = l1.next;
                }
                // similarly if l2 is not null, its value is added to sum
                if (l2 != null) {
                    sum += l2.val;
                    l2 = l2.next;
                }
                sum += carry; // if carry is there, carry is also added to sum
                carry = sum / 10; // carry is calculated, except the last digit, all others of sum are stored as
                                  // carry
                ListNode node = new ListNode(sum % 10); // the last digit of sum is stored as a node in linked list
                temp.next = node; // the temp node will point to the new added node
                temp = temp.next; // temp is made the new node, i.e temp is moved one place ahead (i.e the new
                                  // node)

            }
            return dummy.next;
        }
    }

    // ^ Delete a node when that node is given
    // Write a function to delete a node in a singly-linked list. You will not be
    // given access to the head of the list, instead you will be given access to the
    // node to be deleted directly.
    // It is guaranteed that the node to be deleted is not a tail node in the list.
    // Input: head = [4,5,1,9], node = 5
    // Output:[4,1,9]
    // Explanation: You are given the second node with value 5, the linked list
    // should become 4 -> 1 -> 9 after calling your function.

    class delete_a_node {
        public void deleteNode(ListNode node) {
            node.val = node.next.val; // copy the next node value to the current node
            node.next = node.next.next; // unlink the next node, whose value we just copied
        }
    }

    // ^ Intersection of 2 linked Lists
    // Given the heads of two singly linked-lists headA and headB, return the node
    // at which the two lists intersect. If the two linked lists have no
    // intersection at all, return null.

    public class intersection_of_2_linked_lists_brute {
        public ListNode getIntersectionNode(ListNode headA, ListNode headB) {
            HashSet<ListNode> s = new HashSet<>();
            // add all nodes of first list to the HashSet
            while (headA != null) {
                s.add(headA);
                headA = headA.next;
            }
            // checking if any element from list2 matches
            while (headB != null) {
                if (s.contains(headB))
                    return headB;
                headB = headB.next;
            }
            return null;
        }
    }

    // better approach is to first find the length of both the linked lists
    // then traverse the head of the longer list n spaces, where n is the difference
    // between the
    // lengths of the 2 lists
    // then keep traversing both the the lists until they reach null
    // and keep comparing their nodes, if they match, that point is the intersection
    // point
    // else return null

    public class intersection_of_2_linked_lists_optimal {
        public ListNode getIntersectionNode(ListNode headA, ListNode headB) {
            if (headA == null || headB == null)
                return null;

            // two pointers a and b pointing at head of both linked lists to traverse
            ListNode a = headA;
            ListNode b = headB;

            // if a and b become same, thats the intersection point, be it null also
            while (a != b) {
                // if a becomes null, it will point to head of other linked List, or else it
                // will move one
                // step ahead
                a = a == null ? headB : a.next;
                // same for b also
                b = b == null ? headA : b.next;

                // if both are null at the same time, they will exit the while loop
                // and null will be returned
            }
            return a;
        }
    }

    // ^ Detect a cycle in a Linked List
    // Given head, the head of a linked list, determine if the linked list has a
    // cycle in it.
    // There is a cycle in a linked list if there is some node in the list that can
    // be reached again by continuously following the next pointer. Internally, pos
    // is used to denote the index of the node that tail's next pointer is connected
    // to. Note that pos is not passed as a parameter.
    // Return true if there is a cycle in the linked list. Otherwise, return false.

    // Brute force approach is to use HashSet, add each Node in the set
    // if any Node matches the previous contained Node in the Set, cycle exists,
    // else not

    // optimal approach is to use tortoise method.

    public class detect_cycle_in_linked_list {
        public boolean hasCycle(ListNode head) {
            // if list is empty or haas only one node
            if (head == null || head.next == null)
                return false;

            // we use tortoise method
            // also known as Floyd's cycle detection method
            ListNode slow = head, fast = head;
            do {
                slow = slow.next;
                fast = fast.next.next;
                if (slow == null || fast == null || fast.next == null)
                    return false;
            } while (fast != slow);

            return true;
        }
    }

    // ? Day 7 Linked Lists and Arrays

    // ^ Rotate a Linked List

    // Given the head of a linked list, rotate the list to the right by k places.
    // Example 1:

    // Input: head = [1,2,3,4,5], k = 2
    // Output: [4,5,1,2,3]

    class rotate_list {
        public ListNode rotateRight(ListNode head, int k) {
            // edge cases
            if (head == null || head.next == null || k == 0)
                return head;

            ListNode curr = head;
            int len = 1;
            // calculate the length of the list
            while (curr.next != null) {
                len++;
                curr = curr.next;
            }

            // go to that node
            curr.next = head; // make it a circular linked list
            k = k % len; // len rotations will make it the original list
            k = len - k; // the length from head, where we will break the link

            while (k-- > 0)
                curr = curr.next;

            // make curr's next as head, and unlink
            head = curr.next;
            curr.next = null;

            return head;
        }
    }

    // ^ Clone a Linked List with random and next pointer

    // A linked list of length n is given such that each node contains an additional
    // random pointer, which could point to any node in the list, or null.
    // Example 1:

    // Input: head = [[7,null],[13,0],[11,4],[10,2],[1,0]]
    // Output: [[7,null],[13,0],[11,4],[10,2],[1,0]]

    // ~ Linked List Definition with Random pointer
    class Node {
        int val;
        Node next;
        Node random;

        public Node(int val) {
            this.val = val;
            this.next = null;
            this.random = null;
        }
    }

    class copy_list_with_random_pointer {
        public Node copyRandomList(Node head) {
            Node iter = head, next;

            // First round: make copy of each node,
            // and link them together side-by-side in a single list.
            while (iter != null) {
                next = iter.next;

                Node copy = new Node(iter.val);
                iter.next = copy;
                copy.next = next;

                iter = next;
            }

            // Second round: assign random pointers for the copy nodes.
            iter = head;
            while (iter != null) {
                if (iter.random != null) {
                    iter.next.random = iter.random.next;
                }
                iter = iter.next.next;
            }

            // Third round: restore the original list, and extract the copy list.
            iter = head;
            Node pseudoHead = new Node(0);
            Node copy, copyIter = pseudoHead;

            while (iter != null) {
                next = iter.next.next;

                // extract the copy
                copy = iter.next;
                copyIter.next = copy;
                copyIter = copy;

                // restore the original list
                iter.next = next;

                iter = next;
            }

            return pseudoHead.next;
        }
    }
}
