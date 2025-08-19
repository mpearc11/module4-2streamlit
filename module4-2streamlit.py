#trying to make subprocesses work
import subprocess
import sys

subprocess.run([f"{sys.executable}", "annotate.py"])

#below is sufficient for running locally
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.resources import get_seq_record
from bokeh.io import show
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align
from io import StringIO

aligner = Align.PairwiseAligner()
aligner.mode = 'local'

uploaded_file = st.file_uploader("",type='fasta',key=1)

if uploaded_file is not None:
    st.success("nanopore sequence file uploaded")
else:
    st.info("please upload your nanopore sequence file")

if st.button('annotate sequence'):
    stringio = StringIO(uploaded_file.getvalue().decode('utf-8'))
    record = SeqIO.read(stringio, 'fasta')
    npseq = str(record.seq)

# get pandas df of annotations
hits = annotate(npseq, is_detailed = True, linear= True)

# get biopython SeqRecord object; this is already done earlier i think
#seq_record = get_seq_record(hits, npseq)

# show plot
show(get_bokeh(hits, linear=True))

@st.fragment()
def alignseqs():
    uploaded_file2 = st.file_uploader("",type='fasta',key=2)

    if uploaded_file2 is not None:
        st.success("original target sequence file uploaded")
    else:
        st.info("please upload your target sequence file")

    if st.button('align sequences'):
        stringio = StringIO(uploaded_file2.getvalue().decode('utf-8'))
        record = SeqIO.read(stringio, 'fasta')
        ogseq = str(record.seq)
        alignment = aligner.align(npseq,ogseq)
        st.write(alignment[0])
alignseqs()
