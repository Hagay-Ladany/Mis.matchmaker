#!/usr/bin/perl
use 5.010;
##########user input to variables##########
use CGI qw(:standard);
$first = param ("sss"); #sequence
$last = param ("num");  #position
$avoid = param ("avoid"); #checkbox

## arrays with the enzyme name, and recognition site which is written in regex format , the lists were available at NEB website , I used excel commands to manipulate the text to be in the right format for regex and perl arrays ##


@enzyme_name_on = ("AclI", "HindIII HindIII-HF®", "SspI SspI-HF®", "MluCI", "PciI", "AgeI AgeI-HF®", "BspMI BfuAI", "SexAI", "MluI MluI-HF®", "BceAI", "HpyCH4I", "HpyCH4III", "BaeI", "BsaXI", "AflIII", "SpeI SpeI-HF®", "BsrI", "BmrI", "BglII", "AfeI", "AluI", "StuI", "ScaI-HF®", "BspDI ClaI", "PI-SceI", "NsiI-HF® NsiI", "AseI", "SwaI", "CspCI", "MfeI MfeI-HF®", "BssS?I", "Nb.BssSI", "BmgBI", "PmlI", "DraIII-HF®", "AleI", "EcoP15I", "PvuII PvuII-HF®", "AlwNI", "BtsIMutI", "NdeI", "FatI", "CviAII", "NlaIII", "MslI", "XcmI", "BstXI", "PflMI", "BccI", "NcoI-HF® NcoI", "BseYI", "FauI", "SmaI", "TspMI XmaI", "LpnPI", "AciI", "SacII", "BsrBI", "HpaII MspI", "ScrFI", "StyD4I", "BsaJI", "BslI", "BtgI", "NciI", "AvrII", "MnlI", "Nb.BbvCI", "Nt.BbvCI", "BbvCI", "SbfI SbfI-HF®", "Bpu10I", "Bsu36I", "EcoNI", "HpyAV", "PspGI", "BstNI", "StyI StyI-HF®", "BcgI", "PvuI PvuI-HF®", "BstUI", "EagI EagI-HF®", "RsrII", "BsiEI", "BsiWI-HF® BsiWI", "BsmBI", "Hpy99I", "MspA1I", "SgrAI", "BfaI", "BspCNI", "XhoI PaeR7I", "EarI", "AcuI", "PstI PstI-HF®", "BpmI", "DdeI", "SfcI", "AflII", "BpuEI", "SmlI", "BsoBI AvaI", "MboII", "BbsI BbsI-HF®", "XmnI", "BsmI", "Nb.BsmI", "EcoRI EcoRI-HF®", "HgaI", "ZraI", "AatII", "Tth111I PflFI", "PshAI", "AhdI", "DrdI", "SacI SacI-HF®", "Eco53kI", "BseRI", "Nt.BstNBI", "MlyI", "PleI", "HinfI", "EcoRV EcoRV-HF®", "DpnI", "Sau3AI MboI DpnII", "BsaBI", "TfiI", "Nb.BsrDI", "BsrDI", "BbvI", "Bts?I", "Nb.BtsI", "BstAPI", "SfaNI", "SphI SphI-HF®", "SrfI", "NmeAIII", "NgoMIV", "NaeI", "BglI", "AsiSI", "BtgZI", "HhaI", "HinP1I", "BssHII", "NotI NotI-HF®", "Fnu4HI", "Cac8I", "MwoI", "NheI NheI-HF®", "BmtI BmtI-HF®", "BspQI SapI", "Nt.BspQI", "BlpI", "ApeKI TseI", "Bsp1286I", "AlwI", "Nt.AlwI", "BamHI BamHI-HF®", "BtsCI", "FokI", "HaeIII", "FseI", "SfiI", "NarI", "PluTI", "KasI", "SfoI", "AscI", "EciI", "BsmFI", "PspOMI", "ApaI", "Sau96I", "NlaIV", "KpnI KpnI-HF®", "Acc65I", "BsaI BsaI-HF® BsaI-HF®v2", "HphI", "BstEII  BstEII-HF®", "AvaII", "BanI", "BaeGI", "BsaHI", "BanII", "RsaI", "CviQI", "BstZ17I-HF®", "BciVI", "SalI SalI-HF®", "Nt.BsmAI", "BcoDI BsmAI", "ApaLI", "BsgI", "AccI", "Hpy166II", "Tsp45I", "HpaI", "PmeI", "HincII", "BsiHKAI", "TspRI", "ApoI ApoI-HF", "NspI", "BsrF?I", "BstYI", "HaeII", "CviKI-1", "EcoO109I", "PpuMI", "I-CeuI", "SnaBI", "I-SceI", "BspHI", "BspEI", "MmeI", "Taq?I", "NruI-HF® NruI", "Hpy188I", "Hpy188III", "XbaI", "BclI BclI-HF", "HpyCH4V", "FspI", "PI-PspI", "MscI", "BsrGI-HF® BsrGI", "MseI", "PacI", "PsiI", "BstBI", "DraI", "PspXI", "BsaWI", "BsaAI", "EaeI");
@enzymes_on = ("(AACGTT)","(AAGCTT)","(AATATT)","(AATT)","(ACATGT)","(ACCGGT)","(ACCTGC)","(ACC(A|T)GGT)","(ACGCGT)","(ACGGC)","(ACGT)","(AC.GT)","(AC....GTA(T|C)C)","(AC.....CTCC)","(AC(G|A)(T|C)GT)","(ACTAGT)","(ACTGG)","(ACTGGG)","(AGATCT)","(AGCGCT)","(AGCT)","(AGGCCT)","(AGTACT)","(ATCGAT)","(ATCTATGTCGGGTGCGGAGAAAGAGGTAAT)","(ATGCAT)","(ATTAAT)","(ATTTAAAT)","(CAA.....GTGG)","(CAATTG)","(CACGAG)","(CACGAG)","(CACGTC)","(CACGTG)","(CAC...GTG)","(CAC....GTG)","(CAGCAG)","(CAGCTG)","(CAG...CTG)","(CAGTG)","(CATATG)","(CATG)","(CATG)","(CATG)","(CA(T|C)....(G|A)TG)","(CCA.........TGG)","(CCA......TGG)","(CCA.....TGG)","(CCATC)","(CCATGG)","(CCCAGC)","(CCCGC)","(CCCGGG)","(CCCGGG)","(CC(G|A|T)G)","(CCGC)","(CCGCGG)","(CCGCTC)","(CCGG)","(CC.GG)","(CC.GG)","(CC..GG)","(CC.......GG)","(CC(G|A)(T|C)GG)","(CC(G|C)GG)","(CCTAGG)","(CCTC)","(CCTCAGC)","(CCTCAGC)","(CCTCAGC)","(CCTGCAGG)","(CCT.AGC)","(CCT.AGG)","(CCT.....AGG)","(CCTTC)","(CC(A|T)GG)","(CC(A|T)GG)","(CC(A|T)(A|T)GG)","(CGA......TGC)","(CGATCG)","(CGCG)","(CGGCCG)","(CGG(A|T)CCG)","(CG(G|A)(T|C)CG)","(CGTACG)","(CGTCTC)","(CG(A|T)CG)","(C(A|C)GC(G|T)G)","(C(G|A)CCGG(T|C)G)","(CTAG)","(CTCAG)","(CTCGAG)","(CTCTTC)","(CTGAAG)","(CTGCAG)","(CTGGAG)","(CT.AG)","(CT(G|A)(T|C)AG)","(CTTAAG)","(CTTGAG)","(CT(T|C)(G|A)AG)","(C(T|C)CG(G|A)G)","(GAAGA)","(GAAGAC)","(GAA....TTC)","(GAATGC)","(GAATGC)","(GAATTC)","(GACGC)","(GACGTC)","(GACGTC)","(GAC...GTC)","(GAC....GTC)","(GAC.....GTC)","(GAC......GTC)","(GAGCTC)","(GAGCTC)","(GAGGAG)","(GAGTC)","(GAGTC)","(GAGTC)","(GA.TC)","(GATATC)","(GATC)","(GATC)","(GAT....ATC)","(GA(A|T)TC)","(GCAATG)","(GCAATG)","(GCAGC)","(GCAGTG)","(GCAGTG)","(GCA.....TGC)","(GCATC)","(GCATGC)","(GCCCGGGC)","(GCCGAG)","(GCCGGC)","(GCCGGC)","(GCC.....GGC)","(GCGATCGC)","(GCGATG)","(GCGC)","(GCGC)","(GCGCGC)","(GCGGCCGC)","(GC.GC)","(GC..GC)","(GC.......GC)","(GCTAGC)","(GCTAGC)","(GCTCTTC)","(GCTCTTC)","(GCT.AGC)","(GC(A|T)GC)","(G(G|A|T)GC(A|C|T)C)","(GGATC)","(GGATC)","(GGATCC)","(GGATG)","(GGATG)","(GGCC)","(GGCCGGCC)","(GGCC.....GGCC)","(GGCGCC)","(GGCGCC)","(GGCGCC)","(GGCGCC)","(GGCGCGCC)","(GGCGGA)","(GGGAC)","(GGGCCC)","(GGGCCC)","(GG.CC)","(GG..CC)","(GGTACC)","(GGTACC)","(GGTCTC)","(GGTGA)","(GGT.ACC)","(GG(A|T)CC)","(GG(T|C)(G|A)CC)","(G(G|T)GC(A|C)C)","(G(G|A)CG(T|C)C)","(G(G|A)GC(T|C)C)","(GTAC)","(GTAC)","(GTATAC)","(GTATCC)","(GTCGAC)","(GTCTC)","(GTCTC)","(GTGCAC)","(GTGCAG)","(GT(A|C)(G|T)AC)","(GT..AC)","(GT(G|C)AC)","(GTTAAC)","(GTTTAAAC)","(GT(T|C)(G|A)AC)","(G(A|T)GC(A|T)C)","(..CA(G|C)TG..)","((G|A)AATT(T|C))","((G|A)CATG(T|C))","((G|A)CCGG(T|C))","((G|A)GATC(T|C))","((G|A)GCGC(T|C))","((G|A)GC(T|C))","((G|A)GG.CC(T|C))","((G|A)GG(A|T)CC(T|C))","(TAACTATAACGGTCCTAAGGTAGCGAA)","(TACGTA)","(TAGGGATAACAGGGTAAT)","(TCATGA)","(TCCGGA)","(TCC(G|A)AC)","(TCGA)","(TCGCGA)","(TC.GA)","(TC..GA)","(TCTAGA)","(TGATCA)","(TGCA)","(TGCGCA)","(TGGCAAACAGCTATTATGGGTATTATGGGT)","(TGGCCA)","(TGTACA)","(TTAA)","(TTAATTAA)","(TTATAA)","(TTCGAA)","(TTTAAA)","((G|C|A)CTCGAG(G|T|C))","((A|T)CCGG(A|T))","((T|C)ACGT(G|A))","((T|C)GGCC(G|A))");

@enzyme_name_off = ("AclI", "HindIII HindIII-HF®", "SspI SspI-HF®", "MluCI", "PciI", "AgeI AgeI-HF®", "BspMI BfuAI", "SexAI", "MluI MluI-HF®", "BceAI", "HpyCH4IV", "HpyCH4III", "BaeI", "BsaXI", "AflIII", "SpeI SpeI-HF®", "BsrI", "BmrI", "BglII", "AfeI", "AluI", "StuI", "ScaI-HF®", "BspDI ClaI", "PI-SceI", "NsiI-HF® NsiI", "AseI", "SwaI", "CspCI", "MfeI MfeI-HF®", "BssS?I", "Nb.BssSI", "BmgBI", "PmlI", "DraIII-HF®", "AleI", "EcoP15I", "PvuII PvuII-HF®", "AlwNI", "BtsIMutI", "NdeI", "FatI", "CviAII", "NlaIII", "MslI", "FspEI", "XcmI", "BstXI", "PflMI", "BccI", "NcoI-HF® NcoI", "BseYI", "FauI", "SmaI", "TspMI XmaI", "Nt.CviPII", "LpnPI", "AciI", "SacII", "BsrBI", "HpaII MspI", "ScrFI", "StyD4I", "BsaJI", "BslI", "BtgI", "NciI", "AvrII", "MnlI", "Nb.BbvCI", "Nt.BbvCI", "BbvCI", "SbfI SbfI-HF®", "Bpu10I", "Bsu36I", "EcoNI", "HpyAV", "PspGI", "BstNI", "StyI StyI-HF®", "BcgI", "PvuI PvuI-HF®", "BstUI", "EagI EagI-HF®", "RsrII", "BsiEI", "BsiWI-HF® BsiWI", "BsmBI", "Hpy99I", "MspA1I", "MspJI", "SgrAI", "BfaI", "BspCNI", "XhoI PaeR7I", "EarI", "AcuI", "PstI PstI-HF®", "BpmI", "DdeI", "SfcI", "AflII", "BpuEI", "SmlI", "BsoBI AvaI", "MboII", "BbsI BbsI-HF®", "XmnI", "BsmI", "Nb.BsmI", "EcoRI EcoRI-HF®", "HgaI", "ZraI", "AatII", "Tth111I PflFI", "PshAI", "AhdI", "DrdI", "SacI SacI-HF®", "Eco53kI", "BseRI", "Nt.BstNBI", "MlyI", "PleI", "HinfI", "EcoRV EcoRV-HF®", "DpnI", "Sau3AI MboI DpnII", "BsaBI", "TfiI", "Nb.BsrDI", "BsrDI", "BbvI", "Bts?I", "Nb.BtsI", "BstAPI", "SfaNI", "SphI SphI-HF®", "SrfI", "NmeAIII", "NgoMIV", "NaeI", "BglI", "AsiSI", "BtgZI", "HhaI", "HinP1I", "BssHII", "NotI NotI-HF®", "Fnu4HI", "Cac8I", "MwoI", "NheI NheI-HF®", "BmtI BmtI-HF®", "BspQI SapI", "Nt.BspQI", "BlpI", "ApeKI TseI", "Bsp1286I", "AlwI", "Nt.AlwI", "BamHI BamHI-HF®", "BtsCI", "FokI", "HaeIII", "FseI", "SfiI", "NarI", "PluTI", "KasI", "SfoI", "AscI", "EciI", "BsmFI", "PspOMI", "ApaI", "Sau96I", "NlaIV", "KpnI KpnI-HF®", "Acc65I", "BsaI BsaI-HF® BsaI-HF®v2", "HphI", "BstEII  BstEII-HF®", "AvaII", "BanI", "BaeGI", "BsaHI", "BanII", "RsaI", "CviQI", "BstZ17I-HF®", "BciVI", "SalI SalI-HF®", "Nt.BsmAI", "BcoDI BsmAI", "ApaLI", "BsgI", "AccI", "Hpy166II", "Tsp45I", "HpaI", "PmeI", "HincII", "BsiHKAI", "TspRI", "ApoI ApoI-HF", "NspI", "BsrF?I", "BstYI", "HaeII", "CviKI-1", "EcoO109I", "PpuMI", "I-CeuI", "SnaBI", "I-SceI", "BspHI", "BspEI", "MmeI", "Taq?I", "NruI-HF® NruI", "Hpy188I", "Hpy188III", "XbaI", "BclI BclI-HF", "HpyCH4V", "FspI", "PI-PspI", "MscI", "BsrGI-HF® BsrGI", "MseI", "PacI", "PsiI", "BstBI", "DraI", "PspXI", "BsaWI", "BsaAI", "EaeI");
@enzymes_off = ("(AACGTT)","(AAGCTT)","(AATATT)","(AATT)","(ACATGT)","(ACCGGT)","(ACCTGC)","(ACC(A|T)GGT)","(ACGCGT)","(ACGGC)","(ACGT)","(AC.GT)","(AC....GTA(T|C)C)","(AC.....CTCC)","(AC(G|A)(T|C)GT)","(ACTAGT)","(ACTGG)","(ACTGGG)","(AGATCT)","(AGCGCT)","(AGCT)","(AGGCCT)","(AGTACT)","(ATCGAT)","(ATCTATGTCGGGTGCGGAGAAAGAGGTAAT)","(ATGCAT)","(ATTAAT)","(ATTTAAAT)","(CAA.....GTGG)","(CAATTG)","(CACGAG)","(CACGAG)","(CACGTC)","(CACGTG)","(CAC...GTG)","(CAC....GTG)","(CAGCAG)","(CAGCTG)","(CAG...CTG)","(CAGTG)","(CATATG)","(CATG)","(CATG)","(CATG)","(CA(T|C)....(G|A)TG)","(CC)","(CCA.........TGG)","(CCA......TGG)","(CCA.....TGG)","(CCATC)","(CCATGG)","(CCCAGC)","(CCCGC)","(CCCGGG)","(CCCGGG)","(CC(G|A|T))","(CC(G|A|T)G)","(CCGC)","(CCGCGG)","(CCGCTC)","(CCGG)","(CC.GG)","(CC.GG)","(CC..GG)","(CC.......GG)","(CC(G|A)(T|C)GG)","(CC(G|C)GG)","(CCTAGG)","(CCTC)","(CCTCAGC)","(CCTCAGC)","(CCTCAGC)","(CCTGCAGG)","(CCT.AGC)","(CCT.AGG)","(CCT.....AGG)","(CCTTC)","(CC(A|T)GG)","(CC(A|T)GG)","(CC(A|T)(A|T)GG)","(CGA......TGC)","(CGATCG)","(CGCG)","(CGGCCG)","(CGG(A|T)CCG)","(CG(G|A)(T|C)CG)","(CGTACG)","(CGTCTC)","(CG(A|T)CG)","(C(A|C)GC(G|T)G)","(C..(G|A))","(C(G|A)CCGG(T|C)G)","(CTAG)","(CTCAG)","(CTCGAG)","(CTCTTC)","(CTGAAG)","(CTGCAG)","(CTGGAG)","(CT.AG)","(CT(G|A)(T|C)AG)","(CTTAAG)","(CTTGAG)","(CT(T|C)(G|A)AG)","(C(T|C)CG(G|A)G)","(GAAGA)","(GAAGAC)","(GAA....TTC)","(GAATGC)","(GAATGC)","(GAATTC)","(GACGC)","(GACGTC)","(GACGTC)","(GAC...GTC)","(GAC....GTC)","(GAC.....GTC)","(GAC......GTC)","(GAGCTC)","(GAGCTC)","(GAGGAG)","(GAGTC)","(GAGTC)","(GAGTC)","(GA.TC)","(GATATC)","(GATC)","(GATC)","(GAT....ATC)","(GA(A|T)TC)","(GCAATG)","(GCAATG)","(GCAGC)","(GCAGTG)","(GCAGTG)","(GCA.....TGC)","(GCATC)","(GCATGC)","(GCCCGGGC)","(GCCGAG)","(GCCGGC)","(GCCGGC)","(GCC.....GGC)","(GCGATCGC)","(GCGATG)","(GCGC)","(GCGC)","(GCGCGC)","(GCGGCCGC)","(GC.GC)","(GC..GC)","(GC.......GC)","(GCTAGC)","(GCTAGC)","(GCTCTTC)","(GCTCTTC)","(GCT.AGC)","(GC(A|T)GC)","(G(G|A|T)GC(A|C|T)C)","(GGATC)","(GGATC)","(GGATCC)","(GGATG)","(GGATG)","(GGCC)","(GGCCGGCC)","(GGCC.....GGCC)","(GGCGCC)","(GGCGCC)","(GGCGCC)","(GGCGCC)","(GGCGCGCC)","(GGCGGA)","(GGGAC)","(GGGCCC)","(GGGCCC)","(GG.CC)","(GG..CC)","(GGTACC)","(GGTACC)","(GGTCTC)","(GGTGA)","(GGT.ACC)","(GG(A|T)CC)","(GG(T|C)(G|A)CC)","(G(G|T)GC(A|C)C)","(G(G|A)CG(T|C)C)","(G(G|A)GC(T|C)C)","(GTAC)","(GTAC)","(GTATAC)","(GTATCC)","(GTCGAC)","(GTCTC)","(GTCTC)","(GTGCAC)","(GTGCAG)","(GT(A|C)(G|T)AC)","(GT..AC)","(GT(G|C)AC)","(GTTAAC)","(GTTTAAAC)","(GT(T|C)(G|A)AC)","(G(A|T)GC(A|T)C)","(..CA(G|C)TG..)","((G|A)AATT(T|C))","((G|A)CATG(T|C))","((G|A)CCGG(T|C))","((G|A)GATC(T|C))","((G|A)GCGC(T|C))","((G|A)GC(T|C))","((G|A)GG.CC(T|C))","((G|A)GG(A|T)CC(T|C))","(TAACTATAACGGTCCTAAGGTAGCGAA)","(TACGTA)","(TAGGGATAACAGGGTAAT)","(TCATGA)","(TCCGGA)","(TCC(G|A)AC)","(TCGA)","(TCGCGA)","(TC.GA)","(TC..GA)","(TCTAGA)","(TGATCA)","(TGCA)","(TGCGCA)","(TGGCAAACAGCTATTATGGGTATTATGGGT)","(TGGCCA)","(TGTACA)","(TTAA)","(TTAATTAA)","(TTATAA)","(TTCGAA)","(TTTAAA)","((G|C|A)CTCGAG(G|T|C))","((A|T)CCGG(A|T))","((T|C)ACGT(G|A))","((T|C)GGCC(G|A))");
  
#checkbox to choose between the enzyme lists (checked on = without 2&3 cutters)#
if ($avoid eq 1) {push @enzymes, @enzymes_on;push @enzyme_name, @enzyme_name_on} else {push @enzymes, @enzymes_off;push @enzyme_name, @enzyme_name_off}; 
#######----User sequence Input manipulations---#####
$user_seq = $first;
$SnP = $last;
##########change to capital letters#########
$user_seq =~ (tr/agtc/AGTC/); 
###error messagaes apear when the position is out of sequence##
@seq_o = ("$user_seq");
$n=length($user_seq);
@oneToend =(1..$n);
if ($last > $n) {push @error, "your SNP point is out of sequnce";} 
if ($last <= 0) {push @error, "your SNP point is out of sequnce";}
if ($user_seq eq "") {push @error_seq, "You did not type any sequence";}
#######----correction to user position input and sequence poisition for changes---#######
$SNP = $SnP-1;
$f_change1 = substr $user_seq,($SNP-1),1;
$f_change2 = substr $user_seq,($SNP-2),1;
$f_change3 = substr $user_seq,($SNP-3),1;
$r_change1 = substr $user_seq,($SNP+1),1;
$r_change2 = substr $user_seq,($SNP+2),1;
$r_change3 = substr $user_seq,($SNP+3),1;
# Creation of the sequence variants by using many "if" commands. not so elegant, but works.   ######
foreach $x(0..18){
if ($f_change1 =~ "T") {
                        $var1 = $user_seq;
                        substr($var1,$SNP-1,1)= "C";
                        push @seq, $var1;
                        $var2 = $user_seq;
                        substr($var2,$SNP-1,1)= "G";
                        push @seq, $var2;
                        $var3 = $user_seq;
                        substr($var3,$SNP-1,1)= "A";
                        push @seq, $var3;    
                       }
elsif ($f_change1 =~ "C") {
                        $var1 = $user_seq;
                        substr($var1,$SNP-1,1)= "T";
                        push @seq, $var1;
                        $var2 = $user_seq;
                        substr($var2,$SNP-1,1)= "G";
                        push @seq, $var2;
                        $var3 = $user_seq;
                        substr($var3,$SNP-1,1)= "A";
                        push @seq, $var3;    
                       }
elsif ($f_change1 =~ "G") {
                        $var1 = $user_seq;
                        substr($var1,$SNP-1,1)= "C";
                        push @seq, $var1;
                        $var2 = $user_seq;
                        substr($var2,$SNP-1,1)= "T";
                        push @seq, $var2;
                        $var3 = $user_seq;
                        substr($var3,$SNP-1,1)= "A";
                        push @seq, $var3;    
                       }  
elsif ($f_change1 =~ "A") {
                        $var1 = $user_seq;
                        substr($var1,$SNP-1,1)= "C";
                        push @seq, "$var1";
                        $var2 = $user_seq;
                        substr($var2,$SNP-1,1)= "G";
                        push @seq, $var2;
                        $var3 = $user_seq;
                        substr($var3,$SNP-1,1)= "T";
                        push @seq, $var3;    
                       } ;
if ($f_change2 =~ "T") {
                        $var4 = $user_seq;
                        substr($var4,$SNP-2,1)= "C";
                        push @seq, $var4;    
                        $var5 = $user_seq;
                        substr($var5,$SNP-2,1)= "G";
                        push @seq, $var5;    
                        $var6 = $user_seq;
                        substr($var6,$SNP-2,1)= "A";
                        push @seq, $var6;    
                       }
elsif ($f_change2 =~ "C") {
                        $var4 = $user_seq;
                        substr($var4,$SNP-2,1)= "T";
                        push @seq, $var4;    
                        $var5 = $user_seq;
                        substr($var5,$SNP-2,1)= "G";
                        push @seq, $var5;    
                        $var6 = $user_seq;
                        substr($var6,$SNP-2,1)= "A";
                        push @seq, $var6;    
                       }
elsif ($f_change2 =~ "G") {
                        $var4 = $user_seq;
                        substr($var4,$SNP-2,1)= "C";
                        push @seq, $var4;    
                        $var5 = $user_seq;
                        substr($var5,$SNP-2,1)= "T";
                        push @seq, $var5;    
                        $var6 = $user_seq;
                        substr($var6,$SNP-2,1)= "A";
                        push @seq, $var6;    
                       }  
elsif ($f_change2 =~ "A") {
                        $var4 = $user_seq;
                        substr($var4,$SNP-2,1)= "C";
                        push @seq, $var4;    
                        $var5 = $user_seq;
                        substr($var5,$SNP-2,1)= "G";
                        push @seq, $var5;    
                        $var6 = $user_seq;
                        substr($var6,$SNP-2,1)= "T";
                        push @seq, $var6;    
                       } ;    
if ($f_change3 =~ "T") {
                        $var7 = $user_seq;
                        substr($var7,$SNP-3,1)= "C";
                        push @seq, $var7;    
                        $var8 = $user_seq;
                        substr($var8,$SNP-3,1)= "G";
                        push @seq, $var8;    
                        $var9 = $user_seq;
                        substr($var9,$SNP-3,1)= "A";
                        push @seq, $var9;    
                       }
elsif ($f_change3 =~ "C") {
                        $var7 = $user_seq;
                        substr($var7,$SNP-3,1)= "T";
                        push @seq, $var7;    
                        $var8 = $user_seq;
                        substr($var8,$SNP-3,1)= "G";
                        push @seq, $var8;    
                        $var9 = $user_seq;
                        substr($var9,$SNP-3,1)= "A";
                        push @seq, $var9;    
                       }
elsif ($f_change3 =~ "G") {
                        $var7 = $user_seq;
                        substr($var7,$SNP-3,1)= "C";
                        push @seq, $var7;    
                        $var8 = $user_seq;
                        substr($var8,$SNP-3,1)= "T";
                        push @seq, $var8;    
                        $var9 = $user_seq;
                        substr($var9,$SNP-3,1)= "A";
                        push @seq, $var9;    
                       }  
elsif ($f_change3 =~ "A") {
                        $var7 = $user_seq;
                        substr($var7,$SNP-3,1)= "C";
                        push @seq, $var7;    
                        $var8 = $user_seq;
                        substr($var8,$SNP-3,1)= "G";
                        push @seq, $var8;    
                        $var9 = $user_seq;
                        substr($var9,$SNP-3,1)= "T";
                        push @seq, $var9;    
                       } ;              
                       
if ($r_change1 =~ "T") {
                        $rvar1 = $user_seq;
                        substr($rvar1,$SNP+1,1)= "C";
                        push @seq, $rvar1;    
                        $rvar2 = $user_seq;
                        substr($rvar2,$SNP+1,1)= "G";
                        push @seq, $rvar2;    
                        $rvar3 = $user_seq;
                        substr($rvar3,$SNP+1,1)= "A";
                        push @seq, $rvar3;    
                       }
elsif ($r_change1 =~ "C") {
                        $rvar1 = $user_seq;
                        substr($rvar1,$SNP+1,1)= "T";
                        push @seq, $rvar1;    
                        $rvar2 = $user_seq;
                        substr($rvar2,$SNP+1,1)= "G";
                        push @seq, $rvar2;    
                        $rvar3 = $user_seq;
                        substr($rvar3,$SNP+1,1)= "A";
                        push @seq, $rvar3;    
                       }
elsif ($r_change1 =~ "G") {
                        $rvar1 = $user_seq;
                        substr($rvar1,$SNP+1,1)= "C";
                        push @seq, $rvar1;    
                        $rvar2 = $user_seq;
                        substr($rvar2,$SNP+1,1)= "T";
                        push @seq, $rvar2;    
                        $rvar3 = $user_seq;
                        substr($rvar3,$SNP+1,1)= "A";
                        push @seq, $rvar3;    
                       }  
elsif ($r_change1 =~ "A") {
                        $rvar1 = $user_seq;
                        substr($rvar1,$SNP+1,1)= "C";
                        push @seq, $rvar1;    
                        $rvar2 = $user_seq;
                        substr($rvar2,$SNP+1,1)= "G";
                        push @seq, $rvar2;    
                        $rvar3 = $user_seq;
                        substr($rvar3,$SNP+1,1)= "T";
                        push @seq, $rvar3;    
                       };
if ($r_change2 =~ "T") {
                        $rvar4 = $user_seq;
                        substr($rvar4,$SNP+2,1)= "C";
                        push @seq, $rvar4;    
                        $rvar5 = $user_seq;
                        substr($rvar5,$SNP+2,1)= "G";
                        push @seq, $rvar5;    
                        $rvar6 = $user_seq;
                        substr($rvar6,$SNP+2,1)= "A";
                        push @seq, $rvar6;    
                       }
elsif ($r_change2 =~ "C") {
                        $rvar4 = $user_seq;
                        substr($rvar4,$SNP+2,1)= "T";
                        push @seq, $rvar4;    
                        $rvar5 = $user_seq;
                        substr($rvar5,$SNP+2,1)= "G";
                        push @seq, $rvar5;    
                        $rvar6 = $user_seq;
                        substr($rvar6,$SNP+2,1)= "A";
                        push @seq, $rvar6;    
                       }
elsif ($r_change2 =~ "G") {
                        $rvar4 = $user_seq;
                        substr($rvar4,$SNP+2,1)= "C";
                        push @seq, $rvar4;    
                        $rvar5 = $user_seq;
                        substr($rvar5,$SNP+2,1)= "T";
                        push @seq, $rvar5;    
                        $rvar6 = $user_seq;
                        substr($rvar6,$SNP+2,1)= "A";
                        push @seq, $rvar6;    
                       }  
elsif ($r_change2 =~ "A") {
                        $rvar4 = $user_seq;
                        substr($rvar4,$SNP+2,1)= "C";
                        push @seq, $rvar4;    
                        $rvar5 = $user_seq;
                        substr($rvar5,$SNP+2,1)= "G";
                        push @seq, $rvar5;    
                        $rvar6 = $user_seq;
                        substr($rvar6,$SNP+2,1)= "T";
                        push @seq, $rvar6;    
                       } ;       
if ($r_change3 =~ "T") {
                        $rvar7 = $user_seq;
                        substr($rvar7,$SNP+3,1)= "C";
                        push @seq, $rvar7;    
                        $rvar8 = $user_seq;
                        substr($rvar8,$SNP+3,1)= "G";
                        push @seq, $rvar8;    
                        $rvar9 = $user_seq;
                        substr($rvar9,$SNP+3,1)= "A";
                        push @seq, $rvar9;    
                       }
elsif ($r_change3 =~ "C") {
                        $rvar7 = $user_seq;
                        substr($rvar7,$SNP+3,1)= "T";
                        push @seq, $rvar7;    
                        $rvar8 = $user_seq;
                        substr($rvar8,$SNP+3,1)= "G";
                        push @seq, $rvar8;    
                        $rvar9 = $user_seq;
                        substr($rvar9,$SNP+3,1)= "A";
                        push @seq, $rvar9;    
                       }
elsif ($r_change3 =~ "G") {
                        $rvar7 = $user_seq;
                        substr($rvar7,$SNP+3,1)= "C";
                        push @seq, $rvar7;    
                        $rvar8 = $user_seq;
                        substr($rvar8,$SNP+3,1)= "T";
                        push @seq, $rvar8;    
                        $rvar9 = $user_seq;
                        substr($rvar9,$SNP+3,1)= "A";
                        push @seq, $rvar9;    
                       }  
elsif ($r_change3 =~ "A") {
                        $rvar7 = $user_seq;
                        substr($rvar7,$SNP+3,1)= "C";
                        push @seq, $rvar7;    
                        $rvar8 = $user_seq;
                        substr($rvar8,$SNP+3,1)= "G";
                        push @seq, $rvar8;    
                        $rvar9 = $user_seq;
                        substr($rvar9,$SNP+3,1)= "T";
                        push @seq, $rvar9;    
                       } last;}

@seq22 = @seq;
foreach $jj(0 .. $#seq22) {$JJ=$jj+1; substr (@seq22[$jj], $SnP-1, 0) = '<font color="red">';}
foreach $jj(0 .. $#seq22) {$JJ=$jj+1; substr (@seq22[$jj], $SnP+18, 0) = '</font>';}



foreach $j(0 .. $#seq22) {$J=$j+1; push @var_out, "$J. @seq22[$j]\n"}; #the output is pushed to array###

## foreach loops to search for the enzyme recognition site and the position in sequence. ##
                          foreach (@enzymes){
                                             my $enz_search = qr{$_}i;
                                           foreach $m(0 .. $#seq){   foreach (@seq_o[$m]=~ /($enz_search)/pg){
                                                                     my $start = $-[0];$start_1 = $start+1;
                                                                     my $end = $+[0];
                                                                     my $hitpos = "$start_1-$end";
                                                                     ;push @pozo, " at $hitpos";
                      last }                     }                  ;}
##this part made to extract the enzyme name and pattern
foreach (@enzymes){foreach $i(0 .. $#seq_o){ if (@seq_o[$i] =~ /$_/) { my $search_for = "$_";
                                                                   my( $index )= grep { $enzymes[$_] eq $search_for } 0..$#enzymes;
                                                                   $name_o = @enzyme_name[$index];
                                                                   $pattern_o = "$_ was found in the sequence";
                                                                   push @output_o, "$name_o $pattern_o";  
								   push @elim, "$name_o";           
                                                                  }   
                                         }
                  } ;
## push the resultes to arrays, elim2 array was made to eliminate findings in the variants resultes ##	  
 foreach $q(0 .. $#output_o) {
                                  push @o_output, "@output_o[$q]@pozo[$q]\n";
				  push @elim2, "@elim[$q]@pozo[$q]";
                              }
			      	      
			      
## the same but for the variants #########
                          foreach (@enzymes){
                                             my $enz_search = qr{$_}i;
                                           foreach $m(0 .. $#seq){   foreach (@seq[$m]=~ /($enz_search)/pg){
                                                                     my $start = $-[0];$start_1 = $start+1;
                                                                     my $end = $+[0];
                                                                     my $hitpos = "$start_1-$end";
                                                                     ;push @pozz, " at $hitpos\n";
                      last }                     }                  ;}

                      
foreach (@enzymes){foreach $i(0 .. $#seq){ if (@seq[$i] =~ /$_/) { my $search_for = "$_";
                                                                   my( $index )= grep { $enzymes[$_] eq $search_for } 0..$#enzymes;
                                                                   $name = @enzyme_name[$index];
                                                                   $pattern = "$_ was found in variant";
                                                                   $seq_num = $i+1; 
                                                                   push @output_v_n, "$name $pattern $seq_num";          
                                                                  }   
                                         }
                  } 
my @indexx = sort { @output_v_n[$a] <=> @output_v_n[$b] } 0 .. $#output_v_n;
@output_v_n = @output_v_n[@indexx];
@pozz = @pozz[@indexx];		  
		  
		  
		  
foreach $k(0 .. $#output_v_n){push @output_var, "@output_v_n[$k]@pozz[$k]";}		  

foreach $i(0..$#elim2){push @new, substr @elim2[$i], 0, 5}; #print"array start: @new\n";         #produce a new array - made to to compare the beginings of the original sequence resultes and the variants#####
foreach $i(0..$#elim2){push @endd, substr @elim2[$i],-2}; #print"array end: @endd\n";            #produce a new array - made to to compare the endings of the original sequence and the variants#####

                                                                                                 #To display only new resitriction sites: this code line will search for the "begining" and "ending"in the output of the original sequence and if found in variants do nothing, if not push will push to a new aray 
foreach $i(0..$#output_var){ if ( grep( /^(?=[@new\[$i\]]).*[@endd\[$i\]]$/, @output_var[$i])){
                                                                                               #print ""
											       } else {push @result, @output_var[$i]} };
if (@o_output[0] eq "") {push @o_output, " no findings!"};	## "no findings" messages ##	   
if (@result[0] eq "") {push @result, " no findings!"};	

#eliminate duplications ( I dont know why they apeared!)###########
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
@user_seq_resultes = uniq(@o_output);
@variants_resultes = uniq(@result);


#creation of html #################################    
print "content-type: text/html\n\n";
print "<HTML>\n",
      "<BODY  bgcolor=\"#fffee8\"> \n",;

print "
<head>
<style>
.button {
    background-color: ffeec2;
    border: none;
    color: black;
    border-radius: 6px;
    padding: 13px 28px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 17px;
    margin: 1px 1px;
    cursor: pointer;
    }
.button:hover {
  border: none;
  background: orange;
    color: white;
  box-shadow: 0px 0px 1px #777;
}
</style>
<button onclick=\"goBack()\" class=\"button\" ><b>Go Back</b></button>
<script>
function goBack() {window.history.back();}
</script>
<h2>Your sequence analysis:</h2>
  <div class=\"column\" style=\"background-color:#fffee8;\">\n",;
   print "<b>Your Input sequence:</b><br><b>3'</b> $user_seq<b><font color=\"red\">@error_seq</font> 5' | Selected position:</b> <font size=\"4\">$last</font> <font color=\"red\">@error</font>
<b>| Nucleotide at position: </b>"; print uc substr($user_seq, $SNP, 1);
   print "<BR><b>Variants:</b><br>\n",; 
   
   
   foreach $k(0..$#var_out)  {if ($k=~/^[0-8]$/) {print "0"}; print "@var_out[$k]<br>\n",;}; 
   print "<b>Restriction enzymes found in your sequence:</b><br>\n",;
      foreach $k(0..$#user_seq_resultes)  {print "@user_seq_resultes[$k]<br>\n",;}; 
   print "<b>Restriction enzymes found only in the variants</b>:</b><br>\n",;    
   foreach $k(0..$#variants_resultes)  {print "@variants_resultes[$k]<br>\n",;};
   
print
      "<\/BODY>\n",
      "<\/HTML>\n";