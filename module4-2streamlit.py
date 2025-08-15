from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
from plannotate.resources import get_seq_record
from bokeh.io import show
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

uploaded_file = st.file_uploader("",type='fasta')

if uploaded_file is not None:
    st.success("nanopore sequence file uploaded")
else:
    st.info("please upload your nanopore sequence file")

if st.button('annotate sequence'):
    stringio = StringIO(uploaded_file.getvalue().decode('utf-8'))
    st.write(stringio)
    print(type(stringio))
    #npseq = uploaded_file.read()
    #st.write(npseq)
    record = SeqIO.read(stringio, 'fasta')
    st.write(record)
    npseq = Seq(str(stringio))

    # get pandas df of annotations
    hits = annotate(npseq, is_detailed = True, linear= True)

# get biopython SeqRecord object
seq_record = get_seq_record(hits, npseq)

# show plot
show(get_bokeh(hits, linear=True))
