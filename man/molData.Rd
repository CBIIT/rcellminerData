\name{molData}
\alias{molData}
\docType{data}
\title{
CellMiner molecular profiling data for the NCI-60 cell lines.
}
\description{
  The molData object, of rcellminer package class MolData, organizes CellMiner's
  molecular profiling data for the NCI-60 cancer cell lines. Also provided are cell
  line annotations, including tissue of origin, detailed histology, and source information.
}
\usage{data(molData)}
\details{
Data matrices for the following types of data are provided. NOTE: Missing values are represented as NA. Pubmed IDs (PMID) are provided.

  \itemize{
    \item cop: Copy Number (Combined_aCGH): The gene copy data set values are derived from log 2 probe intensities averaged over 4 platforms: Agilent Human Genome CGH Microarray 44A, Nimblegen HG19 CGH 385K WG Tiling v2.0, Affymetrix GeneChip Human Mapping 500k Array Set, and Illumina Human1Mv1_C Beadchip. For a given gene and cell line, the indicated value is the deviation from the median average intensity value (with the median taken over all genes in that cell line). Values close to zero thus correspond to a copy number of 2, while positive and negative deviations from zero indicate copy gains and losses, respectively. Citation: PMID: 24670534.
    \item mut: Binary Gene Level Summarized Mutation: The binary gene mutation data are derived from the variant level data in exo using the rcellminer function getBinaryMutationData (see data generation script make_rcellminerdata.R for details). The binary gene mutation profiles indicate, for each gene and NCI-60 cell line, the presence or absence of a deleterious (and presumably somatic) mutation. Deleterious mutations are frameshift, nonsense, splice site, or missense mutations predicted to alter function by the SIFT and POLYPHEN2 algorithms. Potential germline mutations are filtered by excluding variants with an estimated frequency greater than 0.005 in the normal population (based the 1000 Genomes and ESP5400 databases).
    \item exo: Exome: Exome regions captured using Agilent SureSelect All Exon v1.0 Kit and sequenced using Illumina Genome Analyzer IIx. Values are percent conversion of genetic variant. This data captures mutations in the NCI-60. The file "exome_gene_probeIds_mapping.txt" is provided, since the CellMiner data dump does not provide a mapping between genes and mutations. Citation: PMID: 23856246. NOTE: 0 indicates no variant allele
    \item exp: Expression: The mRNA transcription level values are average intensity values combined from 4 platforms: Affy HG-U133(A-B), Affy HuEx 1.0, Affy HG-U133 Plus 2.0, Affy HG-U95(A-E). Average intensity values for given gene are presented if all probe expression patterns are well correlated. Citation: PMID: 22802077.
    \item xai: Expression (Avg. log2 Intensity).
    \item pro: Protein (RPPA; Protein Lysate): Selected values are derived from 176 antibodies associated with 94 proteins. For some proteins, expression data is provided from multiple antibodies.  Averaging the latter is not appropriate, but proteins with a single antibody can be selected by excluding probe names with an underscore (yielding a 74 protein data set).  Values are provided as dose interpolation values using a method described in Nishizuka et al., 2003. PMID: 14623978
    \item mir: MicroRNA: 15,000 probes for 723 human and 76 human viral miRNA's. Each slide contains 8 arrays. PMID:20442302.
    \item mda: Metadata (Phenotype Data). Further information: http://discover.nci.nih.gov/cellminer/celllineMetadata.do
      \itemize{
        \item age: Patient Age (Years)
        \item is_epithelial: Is cell derived from epithelial source?
        \item is_p53_mut: Is TP53 mutated?
        \item mdr: Multi-Drug Resistance Assay. PMID: 7969041
        \item doublingtime: Cell doubling time in hours
      }
  }
}
\source{
  CellMiner Website:
  URL: http://discover.nci.nih.gov/cellminer
  Data Set Information:
  URL: http://discover.nci.nih.gov/cellminer/datasetMetadata.do
}
\value{
  the molData object that contains Cellminer molecular profiling information as a MolData object
}
\references{
  Reinhold, W.C., et al. (2012) CellMiner: a web-based suite of genomic and
  pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell
  line set, Cancer research, 72, 3499-3511
}
\concept{
  rcellminerData
}
\keyword{
  datasets
}
