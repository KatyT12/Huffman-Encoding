#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# huffman.py
# (c) Frank Stajano 2022-12-14
# $Id: huffman_template.py 70 2023-02-03 15:50:45Z fms27 $


"""Huffman compression of a sequence of fixed-size symbols.

In this toy implementation, the fixed-size symbols to be encoded are
bytes and the variable-length codewords are bitstrings from the
community-contributed bitstring module (not in Python's standard
library):

available at
  https://github.com/scott-griffiths/bitstring/

documented at
  https://bitstring.readthedocs.io/en/stable/index.html

installable through
  pip install bitstring

Throughout this code we must make and maintain a clear distinction
between bytes and small integers. There is obviously a 1-to-1 mapping
between the integers 0..255 and the possible values for a byte, but a
byte is not an integer and a small integer is not a byte. Our
convention is that we represent the symbols of the sequence to be
encoded as Python bytes at all times. If we sometimes used them as
integers, we would be able to store some lookup tables more compactly
(e.g. a table with the frequencies of each symbol could be just a
256-item list). But the memory gain is insignificant and there are
greater opportunities for bugs when symbols are represented in two
possible ways in the code. So instead we are slightly more wasteful
and use a dictionary (with floats indexed by bytes) to represent that
kind of table.

Note this Python gotcha: it doesn't help that "for s in bytesequence",
where isinstance(bytesequence, bytes), gives int values to s, instead
of the length-1 bytes value that a sensible person would expect (I
iterate over a sequence of BYTES and I get integers instead? Come
on...). Apparently PEP 467 proposes to fix this with "for s in
bytesequence.iterbytes()", where s would indeed take on a length-1
bytes value. But it had not been adopted by mainstream Python at the
time of writing this.
"""


# pylint: disable=invalid-name, misplaced-comparison-constant

import heapq
import bitstring  # see this class's docstring for where to get this


class HuffmanCode:

    """A HuffmanCode object is built, using the well-known greedy
    algorithm, out of a table of expected frequencies for the symbols
    we expect to deal with. Once the object is built, we may invoke
    its encode method on a sequence of bytes, returning an encoded
    (hopefully shorter) sequence of bytes; and we may invoke its
    decode method on a sequence of bytes previously encoded with it,
    recovering the original sequence.

    The HuffmanCode class includes helper utilities (as static
    methods) for adding and removing padding bits (to reversibly and
    economically transform a bitstring of arbitrary length into one
    that fits into an integral number of bytes) and for creating a
    table of occurrences or frequencies given a sequence of symbols.
    """

    def __init__(self, frequencyTable):
        """Take a frequency table (a 256-item dictionary of floats (that add
        up to 1), the floats indexed by all possible byte values,
        giving the expected relative frequencies of each byte in the
        inputs to be encoded. Generate the corresponding Huffman code
        as a PrefixTree (q.v.) and store it as self.tree, an internal
        data structure of this object.

        Note to the candidate: to make a frequency table from a string
        of symbols, you may use the makeOccurrencesTable and
        occurrences2frequencies methods in this class.
        """
        heap = []
        for b in frequencyTable.keys():
            f = frequencyTable[b]
            t = PrefixTree(f,b)
            heap.append(t)
        # Turn into min heap
        m = MinHeapWrapper(heap)
        while m.size() > 1:
            if m .size() == 4:
                print()
            t1 = m.popMin()
            t2 = m.popMin()
            t3 = PrefixTree.fromTwoTrees(t1,t2)
            m.push(t3)
        self.tree = m.popMin()

        symbols = [ bitstring.BitStream()] * 256
        e = bitstring.BitStream()
        code_array = self.search_tree(self.tree,e,symbols)
        self.code_array = code_array
        self.EPSILON = 0.000001

    def search_tree(self,tree,bitseq,symbols):
        if tree.isSingleton():
            s = tree.symbol
            symbols[int.from_bytes(s,"little")] = bitseq
            return symbols
        else:
            # Trees should always be full
            t1 = tree.lChild
            t2 = tree.rChild

            b1 = bitseq.copy()
            b1.append([False])
            b2 = bitseq.copy()
            b2.append([True])
            symbols = self.search_tree(t1,b1, symbols)
            symbols = self.search_tree(t2,b2,symbols)
            return symbols

    def encode(self, plaintextBytes):
        """Take a bytes object (immutable array of bytes) to be
        encoded. Encode it by replacing each of its bytes with the
        corresponding bitstring codeword, according to this object's
        Huffman code. Return a bitstring consisting of the
        concatenation of all these codewords, finished off with
        suitable padding (cfr paddingSuitableFor and removePadding
        methods).
        """
        ret = bitstring.BitStream()
        for byte in plaintextBytes:
            byte = int.to_bytes(byte,1,"little")
            ret.append(self.codewordFor(byte))
        padding = self.paddingSuitableFor(ret)
        ret.append(padding)

        #print(ret.bytes())
        return ret

    def traverseTree(self,codeWord):
        t = self.tree
        while not (t.isSingleton()):
            if not codeWord[0]:
                t = t.lChild
            else:
                t = t.rChild
            codeWord = codeWord[1:]
        return t.symbol,codeWord



    def decode(self, encodedAndPaddedBits):
        """Take a bitstring to be decoded, consisting logically of a sequence
        of codewords followed by padding, but practically of an
        immutable bitstring, i.e. a sequence of bits without any
        delimiters to separate the variable-length codewords. Remove
        the padding and substitute each codeword with its
        corresponding byte symbol, according to this object's Huffman
        code. Return a bytearray with the decoded result.
        """
        encodedBits = self.removePadding(encodedAndPaddedBits)
        byteArray = b''
        while len(encodedBits) > 0:
            symbol,encodedBits = self.traverseTree(encodedBits)
            byteArray += symbol
        return byteArray

    def codewordFor(self, symbol):
        """Given a symbol (a single byte), return a bitstring with the
        codeword for that symbol.
        """
        return self.code_array[int.from_bytes(symbol,"little")]

    @staticmethod
    def paddingSuitableFor(bits : bitstring.BitStream):
        """Take a bitstring of arbitrary length and return a short bitstring
        of padding, of length between 1 and 8 bits that, if appended
        to it, will make the total length a multiple of 8. This
        padding is reversible. It consists of a 1 and then as many 0s
        as necessary to reach the next multiple of 8.
        """
        n = bits.length
        n = 8 - (n - (n // 8)*8)
        ret = bitstring.BitStream()
        if n > 0:
            ret.append([True])
            n -= 1
            while n > 0:
                ret.append([False])
                n -= 1
        return ret


    @staticmethod
    def removePadding(bits):
        """Take a padded bitstring, whose length will be a multiple of
        8. Return a new (mutable) bitstring.BitStream obtained from
        the previous one by removing the padding (without changing the
        original). Take away all consecutive trailing 0s, if any, and
        then the first 1. The returned result will be 1 to 8 bits
        shorter than the input.
        """
        while not bits[-1]:
            bits = bits[:-1]
        bits = bits[:-1]
        return bits

    @staticmethod
    def makeOccurrencesTable(symbols):
        """Take a sequence of symbols (a bytes object in Python, i.e. an
        immutable array of bytes) and return a 256-item dictionary of
        non-negative integers, indexed by symbols (bytes), giving the
        number of occurrences of each byte in the given symbol
        sequence. NB: the dictionary will contain entries for all 256
        possible symbols (bytes), even if they don't all occur in the
        given sequence.
        """
        ret = {}
        for s in symbols:
            s = int.to_bytes(s,1,"little")
            if s in ret.keys():
                v = ret[s]
                ret[s] += 1
            else:
                ret[s] = 1

        for i in range(256):
            b = int.to_bytes(i,1,"little")
            if not (b in ret.keys()):
                ret[b] = 0
        return ret

    @staticmethod
    def occurrences2frequencies(occurrences):
        """Take a table of occurrences, as generated by the
        makeOccurrencesTable method. Return the corresponding table of
        frequencies (a 256-item dictionary of floats indexed by
        symbols) obtained by normalising the entries of the previous
        table so that they all add up to 1.0.

        Raise a ValueError exception if all the occurrences in the
        table were 0, because this makes normalisation impossible (in
        the sense that it makes it impossible for this routine to
        guarantee the "sum is 1" postcondition, so we refuse to
        operate on such degenerate tables).
        """
        total = 0
        for o in occurrences.keys():
            total += occurrences[o]

        print("Hello")
        if total == 0:
            raise ValueError()
        ret = {}
        for o in occurrences:
            v = occurrences[o]/total
            ret[o] = v

        print(type(ret))
        return ret

class PrefixTree:
    def __init__(self,frequency,symbol,lChild =None, rChild = None):
        self.symbol = symbol
        self.frequency = frequency
        self.lChild = lChild
        self.rChild = rChild

    @staticmethod
    def fromTwoTrees(t1,t2):
        root_key = t1.key() + t2.key()
        newTree = PrefixTree(root_key, "", t1,t2)
        return newTree

    def find_smallest(self,smallest):
        if self.isSingleton():
            if int.from_bytes(self.symbol,"little") < smallest:
                return self.symbol
            else:
                return smallest
        else:
            s1 = self.lChild.find_smallest(smallest)
            s2 = self.rChild.find_smallest(smallest)
            smallest = min(s1,s2)
            return smallest

        while(t.isSingleton() == False):
            t = t.lChild
        return t.symbol


    def __lt__(self,other):
        if self.key() == other.key():
            s1 = self.find_smallest(1000)
            s2 = other.find_smallest(1000)
            return s1 < s2
        else:
            return self.key() < other.key()

    def key(self):
        return self.frequency

    def isSingleton(self):
        if self.symbol == "":
            return False
        else:
            return True

class MinHeapWrapper:
    def __init__(self,unsortedTrees = None):
        self.heap = []
        if unsortedTrees != None:
            self.heap = unsortedTrees
        heapq.heapify(self.heap)

    def isEmpty(self):
        return len(self.heap) == 0

    def push(self,t):
        heapq.heappush(self.heap,t)

    def popMin(self):
        return heapq.heappop(self.heap)

    def size(self):
        return len(self.heap)



#frame = HuffmanCode.makeOccurrencesTable(["a".encode(),"b".encode(),"c".encode()])
#print(HuffmanCode.occurrences2frequencies(frame))
"""

table = HuffmanCode.makeOccurrencesTable([int.to_bytes(65,1,"little"),int.to_bytes(66,1,"little"),int.to_bytes(67,1,"little")])
table = HuffmanCode.occurrences2frequencies(table)
print(len(table.keys()))
print(table)
h = HuffmanCode(table)
print(h.encode(b'AAABBBCCCABAB'))
"""