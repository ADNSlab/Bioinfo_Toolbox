set.seed(1234)

DE_analysis(geo_accession = "GSE36083", # find this from the GEO accession page
            gpl_id = "GPL570", # find this from the GEO accession page [Platform]
            gsms = "010101" # determine from the experiment summary + overall design. Consider 0 for Control, 1 for case
            )
