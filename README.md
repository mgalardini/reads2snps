reads2snps
==========

From second generation sequencing reads to variant calling.
Also variant calling using pairwise whole genome alignment.

Usage
-----

To map reads to a target genome and call variants type `make map`.
To run a whole genome alignment (with a mask on regions with low coverage),
type `make align`; this assumes that reads from which de novo assembly has been
derived are present. If they are not present, but a mask is known, just add a
`mask.bed` file to the top-level directory. If neither masks or reads are
available, the whole genome alignment can be run using the `make alignnoreads`.

Prerequisites
-------------

* Reference genome in Fasta and GenBank format
* Map and aligning using reads:
    * fastqc
    * seq_crumbs
    * python and biopython
    * seqtk
    * bwa
    * samtools
    * picard
    * gatk
    * freebayes
    * vcflib
* Whole genome alignment:
    * mummer
    * bedtools
    * parsnp

Copyright
---------

Copyright (C) <2015> EMBL-European Bioinformatics Institute

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
GNU General Public License for more details.

Neither the institution name nor the name reads2snps
can be used to endorse or promote products derived from
this software without prior written permission.
For written permission, please contact <marco@ebi.ac.uk>.

Products derived from this software may not be called reads2snps
nor may reads2snps appear in their names without prior written
permission of the developers. You should have received a copy
of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.
