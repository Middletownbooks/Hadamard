**Follow-up: a Hadamard matrix of order 2092 = 4·523**

Posted 2026-10-14, as a timestamped follow-up to my H(2292) post. As far as I can tell, no Hadamard matrix of order 2092 was known before this. I searched recent arXiv (math.CO) and found no paper constructing one; corrections are welcome.

**What it is.** A Hadamard matrix is a ±1 matrix H with H·Hᵀ = n·I (its rows are mutually orthogonal). It is conjectured that one exists for every n divisible by 4. The smallest orders with no known construction are in the low thousands, and 2092 was one of them.

**Construction.** Below are four ±1 sequences A, B, C, D of length 522. Their periodic autocorrelations sum to −4 at every nonzero shift, and their sums are (2, 0, 0, 0). Put them into the Goethals–Seidel array as circulants (back-diagonal R), then add a 4-row, 4-column border:
- border rows S ⊗ 1ᵀ and border columns Tᵀ ⊗ 1, with S, T ∈ H(4);
- the 4×4 corner −S·Tᵀ/2.

This gives a 2092 × 2092 Hadamard matrix. The S, T and corner are the same as for H(2292).

Every sequence is invariant under multiplication of the index by 181 mod 522. On Z₅₂₂ ≅ Z₁₈ × Z₂₉ this means each sequence is constant on the quartic-residue classes mod 29. That symmetry reduced the search space enough for a two-stage search. The final search ran for about 2½ days on a 4-core machine.

**Credit.** The method was designed, and the search written and run, by Claude Opus 5.5 instances (Anthropic) working in my repository. My role was directing the work and checking it. The result was verified by three separate scripts, including one sharing no code with the search.

    A=---+-++--+-+++-+-+++++--------+--++---+---++-+++-+-+-+---++--+++--+++-----++++-++++-+++--+-+++--+-+--+++-+++-++---++-+++-++---++-+---+--++-+---+--++-++++--+-++--+-----+--+----+--+----++++++-+++---++-------+++--++++++-+------+++-++++-+-+++--++-----+---+-++---+-+-++-+++-+-----++-+-+--+---+++---+++++--+-++++-+-+-++-++-+----+--++-----+----+--+---+-+++-+----+++-+++---++++---+++++--+-+-+--+-++--+-++---+-----+----++++++-+++--+--++-+--++-+--+-+++++-+-++++------+++++++--++-++-++++-++++-+--+--+-+-+++-+--+-+++--+++---+++-++-+--
    B=+-+++---+----+-----+-+++--+-++----+--++-+--+++-++----+-+--++++--+--+-++--+-+-+-+-++--+++----+-----+-+-----++---+--++++--+++++--+--+--++++-+-+++++--++--+--++++++++-++++------+---++-+----++++-+---++++----+--++-++-++-++++-++-+-+-+--+---+-++-++++++-++++++--+++-++--++-+++--++++-+++-+---++++-+---+++-+-+-+----++-+-+-+---+-++-+++-+---+-+--+--++-++-------+-+++++---+---+--+++-++-----++-++-+------++--+++--+--+--+---+-+-++---+--+--+----++-+-++--+++--++++--+---+-++++-+++-++++-+-++--+++--+---+--+--+-----+-++++-++-+++-++-++-+--+-+-
    C=+++-+++--++++--+---+--+++--+-+-+-+-----++-+---+-++-+-+--+--+--++--+-+-+-----+-+-+---+-++--------++-+-+++-+++---+---+--++--+-+--++-+-+----+-+--+++---++++++--+-+-++--+++++-++++++----+++-+--+-+-+++---+-++++++--+-+--+++++++++----+----++-+-+--++-+--+++---+--+-+-+-+++-++-++-++-+++----++-+----+---++-+++-++++++++-++---+-+-+-----+-+--++--++-+-+---+--+-++++----+-++++---++----+----+++-+-++-+++++++--++-++-++--+++-+--++++-+--+-++-++-++-++--+-------+++-----+---+++++-++---+-+-+-+--+++---++--++---++-++-+-+-++---+-+----++-+--++-+-++-
    D=+-++-++++++++-+-+--+++---+++----+------++++-+-+++++-++----+-+-++--+--++-+-+++---+---+++++-++---+++-+-++++--++++--+----++--+++--+--++--+--+-+-+++--+++++--+-++++++++-+-++-+-+++++++------+----+-++--+-+---+-----+-+-++--+-+--++++-+---+---+++------+-++++-++--+-----+++-++-+-++--++---+-++-+--+++-+--+-+-+-+++----+-+-+--++--++---+---+++-++++-+-++++----++++++---+-+--+--+++---------+-----++--+-+++++---+-++-+----++++-+-+--+-+-+++++--++--++++---+--+-++-+-++-----+++-+----++-++---+--++--++++---++---+-+++---++-----++--++-++--++-+--+-

    S = [[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]]
    T = [[1,1,1,-1],[1,1,-1,1],[1,-1,1,1],[-1,1,1,1]]
    corner = [[-1,-1,-1,-1],[-1,-1,1,1],[-1,1,-1,1],[1,-1,-1,1]]

**Quick check (Python + numpy):** for each sequence x (as ±1), sum over the four of `np.dot(x, np.roll(x, -s))` should be 2088 at s = 0 and −4 for every s from 1 to 521.

SHA-256 of the sequence file (`sequences_Z522.txt`, the A–D lines plus a 2-line header): 7e6a2f19149b4163fad79d4f1f257c8c04313bd70692d796cff145c240c2993e
