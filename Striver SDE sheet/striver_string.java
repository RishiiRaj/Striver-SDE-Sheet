//! Striver's SDE sheet all Strings solutions

import java.util.*;

public class striver_string {

    // ? Day 15 Strings

    // ^ Reverse words in a String.

    // A word is defined as a sequence of non-space characters. The words in s will
    // be separated by at least one space.

    // Return a string of the words in reverse order concatenated by a single space.

    // Note that s may contain leading or trailing spaces or multiple spaces between
    // two words. The returned string should only have a single space separating the
    // words. Do not include any extra spaces.

    // Example 1:

    // Input: s = "the sky is blue"
    // Output: "blue is sky the"

    class reverse_words {
        public String reverseWords(String s) {
            StringBuilder sb = new StringBuilder();

            // created an array to store each word
            String[] newArray = s.split(" ");

            // run a reverse for loop
            for (int i = newArray.length - 1; i >= 0; i--) {
                if (!newArray[i].isEmpty())// if new array is not empty
                {
                    sb.append(newArray[i]);// append each word
                    sb.append(" ");// and append white spaces
                }
            }
            return sb.toString().trim(); // return the string trimmed and converted to string
        }
    }

    // ^ longest Palindrome in a String

    // Given a string s, return the longest palindromic substring in s.
    // Example 1:

    // Input: s = "babad"
    // Output: "bab"
    // Explanation: "aba" is also a valid answer.

    // ^Roman to Integer

    // Roman numerals are represented by seven different symbols: I, V, X, L, C, D
    // and M.
    // Example 2:

    // Input: s = "LVIII"
    // Output: 58
    // Explanation: L = 50, V= 5, III = 3.

    class roman_to_integer {
        public int romanToInt(String S) {
            int ans = 0, num = 0;

            // we use right to left approach
            for (int i = S.length() - 1; i >= 0; i--) {
                switch (S.charAt(i)) {
                    case 'I':
                        num = 1;
                        break;
                    case 'V':
                        num = 5;
                        break;
                    case 'X':
                        num = 10;
                        break;
                    case 'L':
                        num = 50;
                        break;
                    case 'C':
                        num = 100;
                        break;
                    case 'D':
                        num = 500;
                        break;
                    case 'M':
                        num = 1000;
                        break;
                }
                /*
                 * So we can avoid the need for an extra variable here. We do run into the case
                 * of repeated numerals causing an issue (ie, "III"), but we can clear that by
                 * multiplying num by any number between 2 and 4 before comparing it to ans,
                 * since the numerals jump in value by increments of at least 5x.
                 */
                if (4 * num < ans)
                    ans -= num;
                else
                    ans += num;
            }
            return ans;
        }
    }

    // ^ Integer to Roman

    // Roman numerals are represented by seven different symbols: I, V, X, L, C, D
    // and M.

    // Example 1:

    // Input: num = 3
    // Output: "III"
    // Explanation: 3 is represented as 3 ones.

    public class integer_to_roman {
        public String intToRoman(int num) {

            int[] values = { 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1 };
            String[] strs = { "M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I" };

            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < values.length; i++) {
                while (num >= values[i]) {
                    num -= values[i];
                    sb.append(strs[i]);
                }
            }
            return sb.toString();
        }
    }

    // ^ Implement ATOI/STRSTR

    // Implement the myAtoi(string s) function, which converts a string to a 32-bit
    // signed integer (similar to C/C++'s atoi function).

    // Only the space character ' ' is considered a whitespace character.
    // Do not ignore any characters other than the leading whitespace or the rest of
    // the string after the digits.

    // Example 1:

    // Input: s = "42"
    // Output: 42
    // Explanation: The underlined characters are what is read in, the caret is the
    // current reader position.
    // Step 1: "42" (no characters read because there is no leading whitespace)
    // Step 2: "42" (no characters read because there is neither a '-' nor '+')
    // Step 3: "42" ("42" is read in)
    // The parsed integer is 42.
    // Since 42 is in the range [-231, 231 - 1], the final result is 42.

    class atoi {
        public int myAtoi(String str) {
            str = str.trim();
            if (str == null || str.length() == 0)
                return 0;//

            char firstChar = str.charAt(0);
            int sign = 1, start = 0, len = str.length();
            long sum = 0;
            if (firstChar == '+') {
                sign = 1;
                start++;
            } else if (firstChar == '-') {
                sign = -1;
                start++;
            }
            for (int i = start; i < len; i++) {
                // if we find first non digit Character, we return the current sum
                if (!Character.isDigit(str.charAt(i)))
                    return (int) sum * sign;
                sum = sum * 10 + str.charAt(i) - '0';

                // take care for overflow case
                if (sign == 1 && sum > Integer.MAX_VALUE)
                    return Integer.MAX_VALUE;
                if (sign == -1 && (-1) * sum < Integer.MIN_VALUE)
                    return Integer.MIN_VALUE;
            }

            return (int) sum * sign;
        }
    }

    // ^ Longest Common Prefix
    // Write a function to find the longest common prefix string amongst an array of
    // strings.
    // If there is no common prefix, return an empty string "".

    // Example 1:

    // Input: strs = ["flower","flow","flight"]
    // Output: "fl"

    class longest_common_prefix {
        public String longestCommonPrefix(String[] strs) {
            if (strs == null || strs.length == 0)
                return "";
            Arrays.sort(strs); // sorts according to length of strings, string with smallest length comes
                               // first

            String first = strs[0]; // find the smallest string in the array
            String last = strs[strs.length - 1]; // find the largest string in the array
            int c = 0;

            while (c < first.length()) { // run a loop till the length of the smallest string
                if (first.charAt(c) == last.charAt(c)) // check for prefix
                    c++;
                else
                    break;
            }
            return (c == 0) ? "" : first.substring(0, c);
        }
    }

    // ^Repeated String Match
    // Given two strings a and b, return the minimum number of times you should
    // repeat string a so that string b is a substring of it. If it is impossible
    // for b​​​​​​ to be a substring of a after repeating it, return -1.

    // Notice: string "abc" repeated 0 times is "", repeated 1 time is "abc" and
    // repeated 2 times is "abcabc".

    // Example 1:

    // Input: a = "abcd", b = "cdabcdab"
    // Output: 3
    // Explanation: We return 3 because by repeating a three times "abcdabcdabcd", b
    // is a substring of it.

    class repeated_string_match {
        public int repeatedStringMatch(String a, String b) {
            StringBuilder sb = new StringBuilder(); // create an object of string builder class
            int count = 0;

            // append a in sb until length of sb exceeds length of b
            while (sb.length() < b.length()) {
                sb.append(a);
                count++;
            }
            if (sb.toString().contains(b))
                return count; // if pattern is found
            if (sb.append(a).toString().contains(b))
                return ++count; // if not, append once again and check
            else
                return -1; // else return -1
        }
    }

    // ^Implement strStr()

    // Given two strings needle and haystack, return the index of the first
    // occurrence of needle in haystack, or -1 if needle is not part of haystack.

    // Example 1:

    // Input: haystack = "hello", needle = "ll"
    // Output: 2

    class str {
        public int strStr(String haystack, String needle) {
            return haystack.indexOf(needle);
        }
    }

    // ^ Check for Anagrams
    // Given two strings s and t, return true if t is an anagram of s, and false
    // otherwise.

    // An Anagram is a word or phrase formed by rearranging the letters of a
    // different word or phrase, typically using all the original letters exactly
    // once.

    // Example 1:

    // Input: s = "anagram", t = "nagaram"
    // Output: true

    class anagram_check {

        // The idea is simple. It creates a size 26 int arrays as buckets for each
        // letter in alphabet. It increments the bucket value with String s and
        // decrement with string t. So if they are anagrams, all buckets should remain
        // with initial value which is zero. So just checking that and return
        public boolean isAnagram(String s, String t) {
            if (s.length() != t.length() || s.length() == 0 || t.length() == 0)
                return false;

            int alphabet[] = new int[26];

            for (int i = 0; i < s.length(); i++) {
                alphabet[s.charAt(i) - 'a']++;
                alphabet[t.charAt(i) - 'a']--;
            }

            for (int i : alphabet) {
                if (i != 0)
                    return false;
            }
            return true;
        }
    }

    // ^ Compare version Numbers
    // Given two version numbers, version1 and version2, compare them.

    // Version numbers consist of one or more revisions joined by a dot '.'. Each
    // revision consists of digits and may contain leading zeros. Every revision
    // contains at least one character. Revisions are 0-indexed from left to right,
    // with the leftmost revision being revision 0, the next revision being revision
    // 1, and so on. For example 2.5.33 and 0.1 are valid version numbers.

    // To compare version numbers, compare their revisions in left-to-right order.
    // Revisions are compared using their integer value ignoring any leading zeros.
    // This means that revisions 1 and 001 are considered equal. If a version number
    // does not specify a revision at an index, then treat the revision as 0. For
    // example, version 1.0 is less than version 1.1 because their revision 0s are
    // the same, but their revision 1s are 0 and 1 respectively, and 0 < 1.

    // Return the following:

    // If version1 < version2, return -1.
    // If version1 > version2, return 1.
    // Otherwise, return 0.

    // Example 1:

    // Input: version1 = "1.01", version2 = "1.001"
    // Output: 0
    // Explanation: Ignoring leading zeroes, both "01" and "001" represent the same
    // integer "1".

    // Example 2:

    // Input: version1 = "1.0", version2 = "1.0.0"
    // Output: 0
    // Explanation: version1 does not specify revision 2, which means it is treated
    // as "0".

    class compare_version_numbers {
        public int compareVersion(String version1, String version2) {

            String[] v1 = version1.split("\\.");
            String[] v2 = version2.split("\\.");

            int max = Math.max(v1.length, v2.length); // finding the max length
            for (int i = 0; i < max; i++) {

                // converting every array element to integer till i is less than that array's
                // length. else making it zero.
                int ver1 = i >= v1.length ? 0 : Integer.parseInt(v1[i]);
                int ver2 = i >= v2.length ? 0 : Integer.parseInt(v2[i]);

                // comparing for each array element
                if (ver1 > ver2)
                    return 1;
                else if (ver2 > ver1)
                    return -1;
            }
            return 0;
        }
    }
}