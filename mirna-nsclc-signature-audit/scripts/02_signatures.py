# Step 2 (imported by later steps): the 13 published miRNA sets and their precursor names.
# Verified from PubMed abstracts retrieved in this session. (name, source, type, miRNAs as reported -> precursor names in the Xena/GDC stem-loop file)
SIGS = [
 ("Li2017",       "Li et al., Oncotarget 2017, 10.18632/oncotarget.15596",      "tissue (TCGA-validated)", ["hsa-mir-21","hsa-mir-30d"]),
 ("Gallach2017",  "Gallach et al., Oncotarget 2017, 10.18632/oncotarget.18603", "tissue",                 ["hsa-mir-21","hsa-mir-188"]),
 ("MacDonagh2021","MacDonagh et al., Transl Lung Cancer Res 2021, 10.21037/tlcr-20-959","tissue",       ["hsa-mir-30a","hsa-mir-30b","hsa-mir-30c-1","hsa-mir-30c-2","hsa-mir-34a","hsa-mir-4286"]),
 ("Hu2010",       "Hu et al., J Clin Oncol 2010, 10.1200/JCO.2009.24.9342",      "serum",                  ["hsa-mir-486-1","hsa-mir-486-2","hsa-mir-30d","hsa-mir-1-1","hsa-mir-1-2","hsa-mir-499"]),
 ("Guo2016",      "Guo et al., Mol Med Rep 2016, 10.3892/mmr.2016.5114",        "serum (advanced NSCLC)", ["hsa-mir-1-1","hsa-mir-1-2","hsa-mir-30d","hsa-mir-221","hsa-mir-486-1","hsa-mir-486-2"]),
 ("Zhang2019",    "Zhang et al., Carcinogenesis 2019, 10.1093/carcin/bgy132",    "serum (advanced NSCLC)", ["hsa-mir-191","hsa-mir-28","hsa-mir-145","hsa-mir-328","hsa-mir-18a"]),
 ("Liu2017",      "Liu et al., Oncotarget 2017, 10.18632/oncotarget.14369",      "plasma exosomal",        ["hsa-mir-23b","hsa-mir-10b","hsa-mir-21"]),
 ("Zhou2020",     "Zhou et al., J Clin Lab Anal 2020, 10.1002/jcla.23505",       "serum",                  ["hsa-mir-942","hsa-mir-601"]),
 ("Xue2020",      "Xue et al., Biomed Res Int 2020, 10.1155/2020/9601876",       "serum",                  ["hsa-mir-1228","hsa-mir-181a-1","hsa-mir-181a-2"]),
 ("Yang2019",     "Yang et al., Int J Biol Sci 2019, 10.7150/ijbs.33986",        "serum",                  ["hsa-mir-146b","hsa-mir-29c"]),
 ("Raponi2009",   "Raponi et al., Cancer Res 2009, 10.1158/0008-5472.CAN-09-0587","tissue (SCC)",          ["hsa-mir-146b"]),
 ("Yu2008",       "Yu et al., Cancer Cell 2008, 10.1016/j.ccr.2007.12.008",     "tissue (RT-PCR)",        ["hsa-let-7a-1","hsa-let-7a-2","hsa-let-7a-3","hsa-mir-221","hsa-mir-137","hsa-mir-372","hsa-mir-182"]),
 ("Effusion2012", "Cell-free malignant effusion study, PLoS One 2012 (PMC3427341)", "pleural effusion",      ["hsa-mir-93","hsa-mir-100","hsa-mir-134","hsa-mir-151","hsa-mir-345"]),
]
