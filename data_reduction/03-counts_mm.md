# From Alignments to a counts table

This document assumes [alignment](./alignment.md) has been completed.

**IF for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy.**

```bash
cp -r /share/workshop/mrnaseq_workshop/jli/rnaseq_example/02-STAR_alignment /share/workshop/mrnaseq_workshop/$USER/rnaseq_example/.
cp  /share/workshop/mrnaseq_workshop/jli/rnaseq_example/summary_star_alignments.txt /share/workshop/mrnaseq_workshop/$USER/rnaseq_example/.
```

In this section, we will collate all of the count data into one file for analysis in R.


### <font color='red'> Start Group Exercise: </font>
- Go through the following commands to generate the collated file for analysis in R.
- Answer the questions at the end of this section.

1. First lets make sure we are where we are supposed to be.

    ```bash
    cd /share/workshop/mrnaseq_workshop/$USER/rnaseq_example
    ```

1. First, go back to your 02-STAR_alignment directory. Let's use a wildcard to list all of the counts files from all of the STAR alignment directories:

    ```bash
    ls -lah 02-STAR_alignment/*/*ReadsPerGene.out.tab
    ```

    Take a look at the beginning of one of these files:

    ```bash
    head 02-STAR_alignment/mouse_110_WT_C/mouse_110_WT_C_ReadsPerGene.out.tab
    ```

    <div class="output"> N_unmapped	215558	215558	215558
    N_multimapping	106596	106596	106596
    N_noFeature	364255	2649569	377560
    N_ambiguous	109260	696	46009
    ENSMUSG00000102693.2	0	0	0
    ENSMUSG00000064842.3	0	0	0
    ENSMUSG00000051951.6	0	0	0
    ENSMUSG00000102851.2	0	0	0
    ENSMUSG00000103377.2	0	0	0
    ENSMUSG00000104017.2	0	0	0
    </div>

    The columns are ID, reads map to either strand, reads mapped to forward strand, and reads mapped to the reverse strand and the first four lines are category totals. In this experiment, it looks like the reads are from the reverse strand, due to the much higher mapping numbers in that column and they similar to reads mapped to either strands. So what we want is just that column of numbers (minus the first four lines), for every one of these files.

1. So let's take one file and figure out how to do that, then we will expand it to all the files. First let's just get the rows we want, i.e. everything but the first four:

    ```bash
    tail -n +5 02-STAR_alignment/mouse_110_WT_C/mouse_110_WT_C_ReadsPerGene.out.tab | head
    ```

    When you give the '-n' option for the 'tail' command a number preceded by a '+' sign, it gives you the entire file starting at the line indicated by the number. In this case, we want to skip the first 4 lines, so we start at line 5. We're piping the command to 'head' just to check that it looks correct. You shouldn't see the first four total lines.

    Now, we want only the fourth column (the counts), and in order to get that we pipe the output of the tail command to the 'cut' command, and then redirect the output to a new file:

    ```bash
    tail -n +5 02-STAR_alignment/mouse_110_WT_C/mouse_110_WT_C_ReadsPerGene.out.tab | cut -f4 > mouse_110_WT_C_ReadsPerGene.out.tab.count
    ```

    Now, mouse_110_WT_C_ReadsPerGene.out.tab.count contains a single column of data... counts for each of the genes for that sample.

1.  Now, we want to do these steps for ALL of the read count files... and to do that we will be using a 'for loop' directly on the command line. First, just run a simple 'for loop' that will print out the names of all the files we want to use:

    ```bash
    for sample in `cat samples.txt`; do echo ${sample}; done
    ```

    This command takes all the files that we listed in step 1 and loops through them, one by one, and for every iteration, assigns the filename to the '${sample}' variable. Also, for every iteration, it runs whatever commands are between the 'do' and 'done'.... and every iteration the value of '${sample}' changes. The semi-colons separate the parts of the loop. The 'echo' command just prints the value of $x to the screen... in this case just the filename. However, instead, we will use our previously created command, but with ${sample} instead of the filename, and adding a few things:

    ```bash
    cd /share/workshop/mrnaseq_workshop/$USER/rnaseq_example
    mkdir 03-Counts
    mkdir 03-Counts/tmp
    for sample in `cat samples.txt`; do \
        echo ${sample}
        cat 02-STAR_alignment/${sample}/${sample}_ReadsPerGene.out.tab | tail -n +5 | cut -f4 > 03-Counts/tmp/${sample}.count
    done
    ```

    After this command, there should be a counts file for every sample, in 03-Counts/tmp.

1. Next, we need to get the columns for the final table. Because all of these files are sorted in the exact same order (by gene ID), we can just use the columns from any of the files:

    ```bash
    tail -n +5 02-STAR_alignment/mouse_110_WT_C/mouse_110_WT_C_ReadsPerGene.out.tab | cut -f1 > 03-Counts/tmp/geneids.txt
    head 03-Counts/tmp/geneids.txt
    ```

    Finally, we want to combine all of these columns together using the 'paste' command, and put it in a temporary file:

    ```bash
    paste 03-Counts/tmp/geneids.txt 03-Counts/tmp/*.count > 03-Counts/tmp/tmp.out
    ```

1. The final step is to create a header of sample names and combine it with the temp file. The header is just all of the sample names separated by tabs. And again, since we pasted the columns in sorted order (wildcards automatically sort in order), the columns just need to be in that same order.

    We take the samples.txt file and pipe that to the sort (to ensure they are in the same order) and then 'paste' command with the '-s' option, which takes a column of values and transposes them into a row, separated by the tab character. And finally, let's put everything together:

    ```bash
    cat <(cat samples.txt | sort | paste -s) 03-Counts/tmp/tmp.out > 03-Counts/rnaseq_workshop_counts.txt
    rm -rf 03-Counts/tmp
    head 03-Counts/rnaseq_workshop_counts.txt
    ```

    <div class="output"> jli@ganesh:/share/workshop/mrnaseq_workshop/jli/rnaseq_example$  head 03-Counts/rnaseq_workshop_counts.txt
    mouse_110_WT_C	mouse_110_WT_NC	mouse_148_WT_C	mouse_148_WT_NC	mouse_158_WT_C	mouse_158_WT_NC	mouse_183_KOMIR150_C	mouse_183_KOMIR150_NC	mouse_198_KOMIR150_C	mouse_198_KOMIR150_NC	mouse_206_KOMIR150_C	mouse_206_KOMIR150_NC	mouse_2670_KOTet3_C	mouse_2670_KOTet3_NC	mouse_7530_KOTet3_C	mouse_7530_KOTet3_NC	mouse_7531_KOTet3_C	mouse_7532_WT_NC	mouse_H510_WT_C	mouse_H510_WT_NC	mouse_H514_WT_C	mouse_H514_WT_NC
    ENSMUSG00000102693.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000064842.3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000051951.6	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000102851.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000103377.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000104017.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000103025.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000089699.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ENSMUSG00000103201.2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    </div>

    And now you have a raw counts file that has a count for every gene, per sample. You will use this file for the next step, which is analysis in R.

1. The alignment exercise and the above exercise to generate the counts table used a subset of the original data. But we are goint to transfer the counts table generated using the full data set and samples.txt to your computer using scp or winSCP, or copy/paste from cat [sometimes doesn't work],  

    In Mac/Linux, users can use scp. Windows users use WinSCP. In a new shell session on my laptop. **NOT logged into tadpole. Replace my [your_username] with your username.**

    ```bash
    mkdir -p ~/rnaseq_workshop
    cd ~/rnaseq_workshop
    scp [your_username]@tadpole.genomecenter.ucdavis.edu:/share/workshop/mrnaseq_workshop/Data/rnaseq_workshop_counts.txt .
    scp [your_username]@tadpole.genomecenter.ucdavis.edu:/share/workshop/mrnaseq_workshop/[your_username]/rnaseq_example/samples.txt .
    ```


**Questions**:
1. Open `rnaseq_workshop_counts.txt` in excel (or excel like application), you may have to move the header column 1 cell to the right, and review as a group. *Anything else worth discussing?*
1. Based on `head 02-STAR_alignment/mouse_110_WT_C/mouse_110_WT_C_ReadsPerGene.out.tab` give some reasoning as to why we are choosing the column we are choosing when creating the counts table.


### <font color='red'> End Group Exercise </font>
