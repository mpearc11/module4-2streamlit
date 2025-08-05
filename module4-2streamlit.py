import base64
import glob
import io
import os
import sys

import numpy as np
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plannotate import __version__ as plannotate_version

from . import resources as rsc
from .annotate import annotate
from .bokeh_plot import get_bokeh


def run_streamlit(args):  # args
    sidebar, cite_fund, images = setup_page()

    inSeq = ""

    uploaded_file = st.file_uploader(
        "Choose a file:", type=rsc.valid_fasta_exts + rsc.valid_genbank_exts
    )

    if uploaded_file is not None:
        name, ext = rsc.get_name_ext(
            uploaded_file.name
        )  # unused name -- could add in

        text_io = io.TextIOWrapper(uploaded_file, encoding="UTF-8")
        text = (
            text_io.read()
        )  # saves this from losing in memory when stream is read
        st.success("File uploaded.")

        inSeq = rsc.validate_file(io.StringIO(text), ext)

if inSeq:
    with st.spinner("Annotating..."):
        linear = st.checkbox("Linear plasmid annotation")
        detailed = st.checkbox("Detailed plasmid annotation")

        with open(rsc.get_resource("templates", "FAQ.html")) as fh:
            faq = fh.read()
        sidebar.markdown(faq + images + cite_fund, unsafe_allow_html=True)

        recordDf = annotate(
            inSeq, args.yaml_file, linear, detailed
        )  # args.blast_db

        if recordDf.empty:
            st.error("No annotations found.")
        else:
            st.markdown("---")
            st.header("Results:")

            st.write(
                "Hover mouse for info, click and drag to pan, scroll wheel to zoom"
            )
            st.bokeh_chart(get_bokeh(recordDf, linear), use_container_width=False)
            if linear:
                st.write(
                    r"\*plasmid is displayed as circular, though pLannotate is treating this as a linear construct"
                )
            if detailed:
                st.write(
                    r"\*\*pLannotate is running in Detailed Annotation mode which can find more hits, though may also find more false positives."
                )

            st.header("Download Annotations:")

            # write and encode gbk for dl
            if option == upload_option and ext in rsc.valid_fasta_exts:
                submitted_fasta = SeqRecord(
                    Seq(inSeq),
                    name=list(SeqIO.parse(io.StringIO(text), "fasta"))[0].id,
                )
                gbk = rsc.get_gbk(recordDf, inSeq, linear, submitted_fasta)

            elif option == upload_option and ext in rsc.valid_genbank_exts:
                submitted_gbk = list(SeqIO.parse(io.StringIO(text), "gb"))[0]
                submitted_gbk.features = []  # clears out old features
                gbk = rsc.get_gbk(recordDf, inSeq, linear, submitted_gbk)

            b64 = base64.b64encode(gbk.encode()).decode()
            gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.gbk"> download {name}_pLann.gbk</a>'
            st.markdown(gbk_dl, unsafe_allow_html=True)

            # encode csv for dl
            cleaned = rsc.get_clean_csv_df(recordDf)
            csv = cleaned.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            csv_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.csv"> download {name}_pLann.csv</a>'
            st.markdown(csv_dl, unsafe_allow_html=True)

            if option == upload_option and ext in rsc.valid_genbank_exts:
                st.header("Download Combined Annotations:")
                st.subheader("uploaded Genbank + pLannotate")
                submitted_gbk = list(SeqIO.parse(io.StringIO(text), "gb"))[0]
                gbk = rsc.get_gbk(recordDf, inSeq, linear, submitted_gbk)
                b64 = base64.b64encode(gbk.encode()).decode()
                gbk_dl = f'<a href="data:text/plain;base64,{b64}" download="{name}_pLann.gbk"> download {name}_pLann.gbk</a>'
                st.markdown(gbk_dl, unsafe_allow_html=True)

            st.markdown("---")

            # prints table of features
            st.header("Features")
            displayColumns = [
                "Feature",
                "percent identity",
                "percent match length",
                "Description",
                "database",
            ]
            markdown = cleaned[displayColumns].copy()
            numericCols = ["percent identity", "percent match length"]
            markdown[numericCols] = np.round(markdown[numericCols], 1)
            markdown[numericCols] = markdown[numericCols].astype(str) + "%"
            markdown.loc[markdown["database"] == "Rfam", "percent identity"] = (
                "-"  # removes percent from Rfam hits
            )
            markdown.loc[markdown["database"] == "Rfam", "percent match length"] = (
                "-"  # removes percent from Rfam hits
            )
            markdown = markdown.set_index("Feature", drop=True)
            markdown = markdown.drop("database", axis=1)
            st.markdown(markdown.drop_duplicates().to_markdown())

            upload_option = 'upload orginal dna fasta file'
            uploaded_file2 = st.file_uploader("upload original dna fasta file",type='fasta')

            if uploaded_file2 is not None:
                nanoseq = SeqRecord(Seq(uploaded_file), id="seq1")
                ogfasta = SeqRecord(Seq(uploaded_file2), id="seq2")
                msa = MultipleSeqAlignment([nanoseq, ogfasta])
                print(msa)
