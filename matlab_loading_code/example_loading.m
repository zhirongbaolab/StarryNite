tstart=1; %start of pharynx
tend=50; %end of retraction


emb_2='L:\examplezipfilelocation\KB_BV460_lron-1_RNAiA_lim-4_ceh-10_08252017_2_w2iSIM-TxRed-600-50_s6_embWT_andmutant_real1.zip'; %e2

templocation='temp_unzip\';
%unzip zipfile to temp file
if ~exist(emb_2,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb_2,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end
%cells_2 is a per cell, cell array of  data
%emb_2  is contents of acetree file organized per timepoint
[ cells_2,emb_e2] = loadcells_unnamed(templocation,tend,1,false );
rmdir(templocation,'s');