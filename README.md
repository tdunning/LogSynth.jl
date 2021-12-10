# LogSynth.jl

A set of random value samplers that help generate fairly realistic
data.

## Skip lists and sampling

A skip list is a data structure in which each entry has different
"heights".  Most entries are low, a few randomly selected entries have
higher heights. Searching in a skip list consists of scanning through
the higher elements to get near, but just before the desired element
and then scanning the next level down from there. This is fast because
we skip over lots of lower elements as we scan the higher ones, hence
the name skip list. This is shown below

<img src="https://user-images.githubusercontent.com/250490/145648549-f2415538-aafe-41f4-be03-689625daf9ea.png" width="400px" style="max-height: 200px; max-width: 100px;"/>

Normally, skip lists are constructed using separately allocated
entries with variable sized arrays of bi-directional links. That makes
scanning efficient in either direction and, importantly, makes
insertion into the middle of the structure fast. A normal skip list
also contains the actual key in this structure.

The skip list distribution data structure is intended to support
querying of a list of values that each have associated weights.  For
our purposes of sampling from a growing set of discrete values, we
don't need to insert new values in the middle of the skip list. Also,
we want to be able to increase the weight of values in the middle of
the list, without having to change a bunch of entries. We also want to
avoid the fragmented memory structure of a normal linked-list
implementation.

To make this work, we turn the skip list on its side and store a
separate vector (really two vectors) for each level. Each entry in
this notional vector includes an index for the same entry at the lower
level and a combined weight for all of the entries before the next
entry at this level. This allows us to append new entries by actually
appending new elements up to the level of the new entry and simply
incrementing the weights of the last entries for all higher
levels. Incrementing the weight of an existing value simply requires
incrementing weights of a single entry for each level. This is shown 
in the following figure.



<img src="https://user-images.githubusercontent.com/250490/145650332-47363e1e-4405-4882-8b93-e6c8c9c88c88.png" width="150px" style="max-height: 200px; max-width: 100px;"/>

We also keep the weights in the skip list un-normalized. To search
with a random value ``u`` in ``[0,1]``, we multiply ``u`` by the total
weight of all samples in the list (which is the sum of the weights in
the highest level vector). We scan at progressively lower levels
keeping a cumulative sum of all the weights we have passed until we
find the particular sample that spans the precise value we want to
find.

The result is O(log n) insertion, update and search times. In
operation, these each take all about 100ns or less on a single laptop
core.

## Sampling from fixed multinomial distributions

If you have a discrete probability distribution with up to a million
discrete values, and where the probability of each value is fixed, the
standard technique is known as an alias table. The most important
features of an alias table are that sampling takes unit time and
memory overhead is very low.

The way that an alias table works is that an alias table with n
entries is created. The entries each have two return values and a
conditional probability. The conditional probability is the
probability of returning the first of the return values given that
this entry is selected.

To create the alias table, the the values to be sampled are sorted
according to probability and two lists are created, one with all
values with probability less than 1/n and the other with all values
greater. An entry is created in the alias table with each value ``i``
that has probability exactly ``1/n``. These entries have both return
values set to ``i`` and a conditional probability set to 1.

Additional entries in this new table is constructed by taking a value
``i`` from the small list and pairing it with a value ``j`` from the
big list.  These new entries have return values ``i``, ``j`` and
conditional probability of ``n p(i)``. The probability of ``j`` is
reduced by ``p(i)``. If the residual probability is still greater than
``1/n``, then ``j`` is left in place and the process repeats with the
next value from the small list. If the residual probability is equal
to ``1/n``, then an entry ``(j, j, 1)`` is added to the alias
table. If the residual probability is less than ``1/n``, then ``j`` is
moved to the small list.

This construction is guaranteed to complete if arithmetic is perfect,
but can be forced to complete if either the big or small lists are
exhausted. Any remaining elements in the other list will have almost
exactly ``1/n`` probability and can simply be treated as if they were
exactly equal.

Sampling consists of generating a random value in
``[0,1)``. An element ``(i, j, p)`` in the alias table is selected
using `ceil(u * n)` and ``i`` is returned if `frac(u * n) < p`.
