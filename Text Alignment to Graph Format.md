#Text Alignment to Graph (TAG) Format Specification

This document defines a file format, the **Text Alignment to Graph** or **TAG** format, to be used by participants in the GA4GH HGVM Graph Aligner Bake-Off. In the bake-off, participants will create Docker containers for their aligners, which take as input a GA4GH graph server URL and a paired-end FASTQ file, and which produce output in TAG format. A TAG file defines **alignments** between a collection of **reads**, grouped into **fragments** (also known as **templates**, generally with two reads each), and a **side graph**. Each read may have zero or more alignments, and all alignments are specified in side graph coordinates.

Here is an example of TAG format.
```
# FRAGNAME	READNUM	FLAGS	PATH	MAPQ	CIGAR	NEXT	SEQ	QUAL
# Two paired end reads from a fragment
frag1	0	99	1:0:5:+,2:0:5:-	30	10=	5:1	GATTACATTA	QQQQQQQQQQ
frag1	1	147	5:1:10:+	30	5=1D3I2=	1:0	GCAGATCCCA	QQQQQQQQQQ
# Two paired unmapped reads
frag2	0	77	*	*	*	*	CAGTAGCCAT	QQQQQQQQQQ
frag2	1	141	*	*	*	*	CGATCACTTA	QQQQQQQQQQ
# Two alignments of the same unpaired read
frag3	*	208	10:0:5:+,10:0:5:+	30	10=	*	TACGTTCACA	QQQQQQQQQQ
frag3	*	464	11:0:6:+,11:100:5:+	0	2=1D2=1X4=	*	TACGTTCACA	QQQQQQQQQQ
# Not shown but possible edge cases:
# Multiple alignments of one read in a fragment, with a single alignment of the other read
# Three or more reads in a fragment
# Fragments with some reads aligned and others unaligned
```

##Overview

TAG is a line-oriented, tab-separated format, in which, as in SAM, `*` is used to represent empty field values. Lines beginning with `#` are comments. Lines beginning with `@` are interpreted as header lines, and are allowed, but their behavior is implementation-defined.

Each line in a TAG file represents an alignment. (Note that some alignments--those with no `PATH` and flagged as unmapped--represent a lack of alignment to anything in the reference, and should be produced for a read if and only if no other alignments are produced. Each alignment is associated with a particular read, and that read is associated with a particular fragment.

##Field definitions

The fields are defined as follows:

* **FRAGNAME**: Name of the fragment that the alignment's read belongs to. Together with `READNUM`, uniquely identifies the read being aligned. May not be `*` or begin with `#`.

* **READNUM**: 0-based index of the read along its fragment. For paired-end reads, must be 0 for read 1 and 1 for read 2. For reads alone on their fragments, may be 0 or `*`.

* **FLAGS**: Flags, represented as a decimal number. The available flags (based off of those in SAM) are:
    * 0x1: the alignment's read's fragment has multiple reads. If set, `READNUM` must not be `*`. If unset, flags `0x2`, `0x8`, `0x20`, `0x40`, and `0x80` are ignored.
    * 0x2: all of the alignments on this alignment's read's fragment are aligned "properly" (i.e. mapped in accordance with the ordering, spacing, and orientation that their reads appear in on the fragment).
    * 0x4: this alignment represents an unmapped read. This alignment must be the only alignment for its read.
    * 0x8: the next read in this alignment's read's fragment is unmapped.
    * 0x10: the `SEQ` for this alignment is the reverse complement of the original input sequence for its read. The `PATH` and `CIGAR` fields are oriented to correspond with the `SEQ` provided here. (Note that this is SAM's reverse-strand flag, but in a potentially cyclic graph reference there is no overall reverse strand.)
    * 0x20: the `SEQ` for the primary alignment of the next read in this alignment's fragment is given on the reverse strand. In other words, the primary alignment for the next read on the fragment has flag `0x10` set.
    * 0x40: there are reads after this alignment's read in its fragment.
    * 0x80: there are reads before this alignment's read in its fragment.
    * 0x100: this alignment is secondary (i.e. a multi-mapping). Each read must have exactly one alignment with this flag un-set.
    * 0x200: this alignment has failed filters or quality control and should not be used.
    * 0x400: this alignment's read is a PCR or optical duplicate.
    * 0x800: Ignored. This is the SAM supplementary alignment flag, but in TAG format a split read can be represented in a single alignment.
    
* **PATH**: Path in the graph to which the read is aligned. May be `*` if the read is not mapped. The path takes the form of a comma-separated list of colon-separated 4-tuples, each of which specifies a subsequence in a sequence in the reference side graph, taken in a particular orientation. Each 4-tuple contains:
    * The ID of a sequence.
    * The 0-based leftmost base along the sequence to include in the subsequence.
    * The number of bases to include in the subsequence.
    * The orientation in which to take the subsequence. `+` for its forward orientation, `-` for its reverse-complement orientation.
    
The subsequences specified are attached together end-to-start in the order listed, to define the linear sequence against which the alignment aligns the read, for the purpose of interpreting the `CIAGR` field. Note that 4-tuples in this list must be maximal: if two adjacent subsequences could be described as a single subsequence, they must be. Also note that the alignment is required to be oriented in the orientation that makes its path numerically smallest (with corresponding fields being compared, and `-` being less than `+`). For example, a path `1:10:5:-,1:0:5:-` must instead instead be specified as `1:0:5:+,1:10:5:+`, and the `SEQ`, `QUAL`, `CIGAR` and `FLAGS` fields must be reverse-complemented, reversed, and adjusted, as appropriate, to agree with the numerically smaller path. Finally, note that `PATH` can express arbitrary jumps between pieces of the reference, eliminating the need for both `D` CIGAR operations and split, supplemental alignments.

* **MAPQ**: Mapping quality score as a decimal value. May be `*` if the read is not mapped, or if mapping qualities are not calculated by the aligner. It is -10 * log10(probability that the `PATH` provided is not the path that produced this read). It should account correctly for the aligner's error model (i.e., the case that, instead of being generated by the given path, the read was generated from a different path and was read with one or more errors) and for multimapping (i.e. the case that, instead of being generated by the given path, the read was generated by a different path that spells out the same sequence).

* **CIGAR**: CIGAR string describing the alignment between `SEQ` and the sequence produced by evaluating `PATH` against the graph. Some operations generally allowed in CIGAR strings are not permitted. The available operations are:
    * **I**: an insetrion in `SEQ` relative to the path.
    * **S**: a "soft clipped" sequence present in `SEQ` but not the path. An `S` operation may only appear adjacent to either the end of the CIGAR string or an `H` operation.
    * **H**: a "hard clipped" sequence not present in `SEQ` or the path. An `H` operation may only appear adjacent to the end of the CIGAR string.
    * **=**: a sequence match between `SEQ` and the path.
    * **X**: a sequence mismatch between `SEQ` and the path.
    
Operations which are explicitly disallowed include `M` (which should be replaced with `=` or `X` as appropriate), and `D`, which should be avoided by skipping the deleted bases in the specification of `PATH` for the alignment. Note also that operations must be maximal: adjacent operations of the same type must be combined.

* **NEXT**: If there is a next read on the alignment's read's template, and if that next read is aligned, gives the smallest-ranked base on the `PATH` of the primary alignment of that next read, as a colon-separated sequence:base tuple. Otherwise, must be `*`. Bases are ranked first by sequence number, and then by offset along the sequence. For example, if the next read's primary alignment's `PATH` is `5:10:5:-,4:100:2:-`, this field would be `4:100`. This field is necessary to efficiently find the primary alignment of the next read, when the file is sorted by position. Note that, for the purposes of this field, the fragment is treated as **circular**, so the "next" read of the last read is the first read.

* **SEQ**: The sequence of the alignment's read, in the orientation specified by `FLAGS`, which is also the correct orientation for interpreting `CIGAR` against the path in `PATH`. Must be expressed in upper case. Must not be `*`.

* **QUAL**: The quality scores for the bases given in `SEQ`, in the same orientation as `SEQ`, expressed in FASTQ-alike phred+33 format. May be `*` if qualities are unavailable.

##Sort order

A TAG file may be unsorted, sorted by position, or sorted by read.

###By Position

If a TAG file is sorted by position, the sort key for each alignment is a 2-tuple, consisting of the sequence ID and offset in the sequence of the base covered by the alignment's `PATH` that would produce the smallest such 2-tuple. For example, an alignment with a `PATH` of `10:3:5:-,2:99:30:+,3:5:15:-` would have a sort key of (2, 99), sorting before (3, 0) or (2, 100) and after (2, 98) or (1, 1000).

Alignments for unmapped reads are moved to the end of the file, and sorted by read, as described below.

###By Read
If a TAG file is sorted by read, the `FRAGNAME` and `READNUM` fields, as a 2-tuple, define the sort order. All the alignments for a read thus sort together. Sorting of the alignments for a given read is implementation-defined.

##Notable Differences from SAM

TAG format is partially based on SAM format, with several important changes:

1. Mapping is now to a side graph, instead of a linear reference. Mapping positions have been replaced by paths in the side graph.
2. Split alignments are now to be expressed on a single line, through the use of the `PATH` field.
3. Deletion operations are no longer allowed in the `CIGAR` strings, and are also to be specified via the `PATH`.
4. The generic `M` operation is no longer allowed in the `CIGAR` field, to be replaced by `=` and `X`.
5. Sort order, when sorting by position, is now by the smallest-valued position visited along the `PATH`, and not by the leftmost position on a linear reference contig. This also 
6. A `READNUM` field has been added, because otherwise the relative order of unmapped reads on a fragment with 4 or more reads is unexpressable.

