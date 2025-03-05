def min_operations(S, T, k1, k2):
    len_S = len(S)
    len_T = len(T)

    # Initialize DP table with infinity
    dp = [[float('inf')] * (len_T + 1) for _ in range(len_S + 1)]

    # Base case: Empty S to empty T requires 0 operations
    dp[0][0] = 0

    # Fill first column: Deleting entire S in chunks of k
    for i in range(k1, len_S + 1):
        for k in range(k1, min(k2, i) + 1):
            dp[i][0] = min(dp[i][0], dp[i - k][0] + 1)

    # Fill the rest of the table
    for i in range(len_S + 1):
        for j in range(1, len_T + 1):
            # Case 1: No operation needed (characters match)
            if i > 0 and j > 0 and S[i - 1] == T[j - 1]:
                dp[i][j] = min(dp[i][j], dp[i - 1][j - 1])

            # Case 2: Tandem Duplication (extend S beyond its original length)
            for k in range(k1, k2 + 1):
                if j >= k:
                    # Ensure we have enough characters to duplicate
                    if i >= k and S[i - k:i] == T[j - k:j]:
                        dp[i][j] = min(dp[i][j], dp[i - k][j - k] + 1)

                    # Allow chaining tandem duplications even if `i == len_S`
                    if i == len_S:
                        prev_idx = j - k  # The segment of T that we are trying to duplicate
                        if prev_idx >= k and T[prev_idx - k:prev_idx] == T[prev_idx: j]:
                            dp[i][j] = min(dp[i][j], dp[i][prev_idx] + 1)

            # Case 3: Deletion (remove a segment of length k from S)
            for k in range(k1, k2 + 1):
                if i >= k:
                    dp[i][j] = min(dp[i][j], dp[i - k][j] + 1)

    print(dp)
    return dp[len_S][len_T] if dp[len_S][len_T] != float('inf') else -1  # Return -1 if impossible


print(min_operations("AB", "ABABAB", 2, 3))