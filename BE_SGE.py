#%%
import pandas as pd
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor

#%%

class General:
    def __init__(self):
        self.codon = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu','CTA': 'Leu',
            'CTG': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'GTT': 'Val',
            'GTC': 'Val','GTA': 'Val', 'GTG': 'Val', 'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
            'CCT': 'Pro','CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr',
            'ACG': 'Thr','GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'TAT': 'Tyr', 'TAC': 'Tyr',
            'TAA': 'STOP','TAG': 'STOP', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln', 'AAT': 'Asn',
            'AAC': 'Asn','AAA': 'Lys', 'AAG': 'Lys', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
            'TGT': 'Cys','TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
            'CGG': 'Arg','AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GGT': 'Gly', 'GGC': 'Gly',
            'GGA': 'Gly', 'GGG': 'Gly'}
        self.codon_short = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 
                        'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V','STOP':'*'}
        self.codon_abb = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
                       'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
                       'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
                       'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
                       'GTT': 'V', 'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W',
                       'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

    def read_fasta(self,acc):
        hg19 = pd.read_csv('EssentialData/hg19_Accession.txt',sep='\t',index_col=1)
        hg38 = pd.read_csv('EssentialData/hg38_Accession.txt',sep='\t',index_col=1)

        if acc in list(hg19.index):
            num = hg19.loc[acc,'Chromosome'].split()[1]        
            fasta = self._fasta19(num)
            return fasta
        elif acc in list(hg38.index):   
            num = hg38.loc[acc,'Chromosome'].split()[1]           
            fasta = self._fasta38(num)
            return fasta
        else:
            return 'FALSE'
        

        
    def _fasta38(self,num):
        total = ''
        with open('EssentialData/hg38/chr{}.fa'.format(num) , 'r') as fasta:
            a = fasta.readlines()

            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1

        return total

    def _fasta19(self,num):
        total = ''
        with open('EssentialData/hg19/hg19_{}.fasta'.format(num) , 'r') as fasta:
            a = fasta.readlines()
            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1
        return total


    def rt_sequence(self,seq):
        rt = ''
        dic = {'A':'T','G':'C','T':'A','C':'G','N':'N','a':'t','g':'c','t':'a','c':'g'}
        for i in seq[::-1]:
            rt += dic[i]
        return rt

class extract_cds_from_all_variants:
    def _extract_all_cds(self,eachdf):
        final = {}
        for i in eachdf.index:
            nc_accession = eachdf.loc[i,'nc_accession']
            ccds_id = eachdf.loc[i,'ccds_id']
            
            fasta = General().read_fasta(nc_accession)
            
            strand = eachdf.loc[i,'cds_strand']

            cds_loc_lst = eachdf.loc[i,'cds_locations']
            cds_loc_lst = cds_loc_lst[1:-1]
            cds_loc_lst = cds_loc_lst.split(', ')

            dic = {}
                
            for idx,cds in enumerate(cds_loc_lst):
                cds = cds.split('-')
                start = int(cds[0])
                end = int(cds[1])                

                if strand == '-':
                    read = fasta[start:end+1]
                    read = General().rt_sequence(read)
                    idx = len(cds_loc_lst)-idx
                    dic[idx] = [start,end]
                    dic[idx].append(read)

                    extended = fasta[start-50:end+1+50]
                    extended = General().rt_sequence(extended)

                    dic[idx].append(extended)

                else:
                    read = fasta[start:end+1]
                    dic[idx+1] = [start,end]
                    dic[idx+1].append(read)

                    extended = fasta[start-50:end+51]
                    dic[idx+1].append(extended)

            final[ccds_id] = dic
        return final

    def _make_each_variant_df(self,eachvar_dic):
        n = 1
        k = 0
        newdic = {}
        while n < len(eachvar_dic)+1:
            cds_info = eachvar_dic[n]
            seq = cds_info[2]
            if len(seq) < 3:
                return pd.DataFrame()
            else:
                if n == 1:
                    if seq[:3] != 'ATG':
                        return pd.DataFrame()
                    else:
                        pass
                else:
                    pass


            l = k
            if l == 0:
                true_cds = seq
            else:
                true_cds = seq[3-l:]

            k = self._check_cds_remnant(seq,k)
            if k == 0:
                true_cds = true_cds
            else:
                true_cds = true_cds[:-k]

            cds_info.append(k)
            cds_info.append(true_cds)
            
            newdic[n] = cds_info
            n += 1

        lst = []
        for i in newdic:
            each_info = eachvar_dic[i]
            lst.append(each_info)
        df = pd.DataFrame(lst,columns = ['start','end','CDS_original','CDS_extended','Remnant_NT','trimmed_CDS'])
        return df

    def _check_cds_remnant(self,seq,prev_rem):
        seq = seq[3-prev_rem:]
        rem = len(seq)%3
        return rem 

    def cds_extractor(self,eachdf):
        p53 = self._extract_all_cds(eachdf)
        dic = {}
        for var in p53:
            eachvar_dic = p53[var]
            eachvardf = self._make_each_variant_df(eachvar_dic)
            dic[var] = eachvardf

        return dic
#%%
class EditingPattern_NGG:
    ## C to T, A to G
    ## G to A, T to C 
    def _check_editable_AorC(self,edit_target_seq):
        lst = []
        for i in edit_target_seq:
            if i == 'A':
                lst.append(i)
            elif i == 'C':
                lst.append(i)
            else:
                pass
        return len(lst)

    def _check_editable_GorT(self,edit_target_seq):
        lst = []
        for i in edit_target_seq:
            if i == 'G':
                lst.append(i)
            elif i == 'T':
                lst.append(i)
            else:
                pass
        return len(lst)
    def _edit_every_AorC(self,a_lst,codon_seq,wt_nt,edited_nt):
        edited_codon = ''
        wt_codon = ''
        lst = []
        if len(a_lst) == 1:
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            

        elif len(a_lst) == 2:
            ## edit all 'A'
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            ## edit one 'A'
            wt_codon = ''
            edited_codon = ''
            for idx, nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()     
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])

        else: ## len(a_lst) == 3
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()  
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower() 
                elif idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                           
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                elif idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                     
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                elif idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                          
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
        
        return lst

    def _edit_all_nt_in_one_codon(self,codon_seq,codon_num):
        a_lst = []
        c_lst = []

        for idx,nt in enumerate(codon_seq):
            if nt == 'A':
                a_lst.append(idx)
            elif nt == 'C':
                c_lst.append(idx)
            else:
                pass
        if len(a_lst) ==0:
            if len(c_lst) == 0:
                return pd.DataFrame()
            else:
                ## C 만 있음
                edited_df = self._edit_every_AorC(c_lst,codon_seq,'C','T')
                edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                edited_df['BaseEditor'] = 'CBE'
        else:
            if len(c_lst) == 0:
                edited_df = self._edit_every_AorC(a_lst,codon_seq,'A','G')
                edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                
                edited_df['BaseEditor'] = 'ABE'
            else: ## 둘 다 있음

                cbe_df = self._edit_every_AorC(c_lst,codon_seq,'C','T')
                cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                cbe_df['BaseEditor'] = 'CBE'

                abe_df = self._edit_every_AorC(a_lst,codon_seq,'A','G')
                abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                abe_df['BaseEditor'] = 'ABE'

                edited_df = pd.concat([cbe_df,abe_df])

        edited_df['Codon_Num'] = codon_num
        return edited_df


    def _edit_partial_nt_in_one_codon(self,codon_seq,codon_num,tup,direction):
        
        final = []
        if direction=='left':
            ## remnant가 왼쪽
            editable_seq_in_codon = codon_seq[tup[0]:]
            remnant =  codon_seq[:tup[0]]
            lst = []
            if len(editable_seq_in_codon) == 1:
                if editable_seq_in_codon == 'A':
                    wt = remnant+'a'
                    edited = remnant+'g'
                    lst.append([wt,edited,'ABE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                elif editable_seq_in_codon == 'C':
                    wt = remnant+'c'
                    edited = remnant+'t'
                    lst.append([wt,edited,'CBE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                else:
                    edited_df = pd.DataFrame()


            elif len(editable_seq_in_codon) == 2:
                a_lst = []
                c_lst = []

                for idx,nt in enumerate(editable_seq_in_codon):
                    if nt == 'A':
                        a_lst.append(idx)
                    elif nt == 'C':
                        c_lst.append(idx)
                    else:
                        pass
                
                if len(a_lst) ==0:
                    if len(c_lst) == 0:
                        edited_df = pd.DataFrame()

                    else:
                        ## C 만 있음
                        edited_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'C','T')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        edited_df['BaseEditor'] = 'CBE'
                else:
                    if len(c_lst) == 0:
                        edited_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'A','G')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        
                        edited_df['BaseEditor'] = 'ABE'
                    else: ## 둘 다 있음
                        cbe_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'C','T')
                        cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                        cbe_df['BaseEditor'] = 'CBE'

                        abe_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'A','G')
                        abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                        abe_df['BaseEditor'] = 'ABE'

                        edited_df = pd.concat([cbe_df,abe_df])
                        
                edited_df.index = [i for i in range(edited_df.shape[0])]
                for idx in range(edited_df.shape[0]):
                    wt = edited_df.loc[idx,'WT']
                    edited = edited_df.loc[idx,'Edited']
                    edited_df.loc[idx,'WT'] = remnant+wt
                    edited_df.loc[idx,'Edited'] = remnant+edited

            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)
                        
        else:
            edited_df = pd.DataFrame()
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df) 
        
 
        if direction=='right':
            ## Remnant r가 오른쪽
            editable_seq_in_codon = codon_seq[:-tup[1]]
            remnant =  codon_seq[-tup[1]:]
            lst = []
            if len(editable_seq_in_codon) == 1:
                if editable_seq_in_codon == 'A':
                    wt = 'a'+remnant
                    edited = 'g'+remnant
                    lst.append([wt,edited,'ABE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                elif editable_seq_in_codon == 'C':
                    wt = 'c'+remnant
                    edited = 't'+remnant
                    lst.append([wt,edited,'CBE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                else:
                    edited_df = pd.DataFrame()


            elif len(editable_seq_in_codon) == 2:
                a_lst = []
                c_lst = []

                for idx,nt in enumerate(editable_seq_in_codon):
                    if nt == 'A':
                        a_lst.append(idx)
                    elif nt == 'C':
                        c_lst.append(idx)
                    else:
                        pass
                
                if len(a_lst) ==0:
                    if len(c_lst) == 0:
                        edited_df = pd.DataFrame()

                    else:
                        ## C 만 있음
                        edited_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'C','T')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        edited_df['BaseEditor'] = 'CBE'
                else:
                    if len(c_lst) == 0:
                        edited_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'A','G')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        
                        edited_df['BaseEditor'] = 'ABE'
                    else: ## 둘 다 있음
                        cbe_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'C','T')
                        cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                        cbe_df['BaseEditor'] = 'CBE'

                        abe_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'A','G')
                        abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                        abe_df['BaseEditor'] = 'ABE'

                        edited_df = pd.concat([cbe_df,abe_df])

                edited_df.index = [i for i in range(edited_df.shape[0])]
                for idx in edited_df.index:
                    wt = edited_df.loc[idx,'WT']
                    edited = edited_df.loc[idx,'Edited']
                    edited_df.loc[idx,'WT'] = wt+remnant
                    edited_df.loc[idx,'Edited'] = edited+remnant
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)

        else:
            edited_df = pd.DataFrame()
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)
        final = pd.concat(final)
        return final

    
    def possible_editing_pattern(self,tup,editable_codon,edit_window):
        ## check 'A' or 'C' in edit_target_seq
        num_editable_codon = len(editable_codon)
        lst = []
        if self._check_editable_AorC(edit_window) == 0:
            return pd.DataFrame()
        else:
            pass
        
        for i,idx in enumerate(editable_codon):
            codon_tup = editable_codon[idx]
            codon_seq = codon_tup[1]
            codon_num = codon_tup[0]

            if self._check_editable_AorC(codon_seq) == 0:
                return pd.DataFrame()
            else:
                pass
            if tup[0] > 0:
                if i == 0:
                    eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'left')

                    lst.append(eachdf)
                else:
                    if tup[1] > 0:
                        
                        if i == num_editable_codon-1:
                            eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'right')
                            lst.append(eachdf)
                        else:
                            eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                            lst.append(eachdf)
                    else:
                        eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                        lst.append(eachdf)
                
            else: ## tup[0] <=0
                if tup[1] > 0:
                    if i == num_editable_codon-1:
                        eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'right')
                        lst.append(eachdf)
                    else:
                        eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                        lst.append(eachdf)
                else:
                    eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                    lst.append(eachdf)

        
        final = pd.concat(lst)
        return final

        

  

        
 

class EditingPattern_CCN:
    
    ## G to A, T to C 
    
    def _check_editable_GorT(self,edit_target_seq):
        lst = []
        for i in edit_target_seq:
            if i == 'G':
                lst.append(i)
            elif i == 'T':
                lst.append(i)
            else:
                pass
        return len(lst)


    def _edit_every_AorC(self,a_lst,codon_seq,wt_nt,edited_nt):
        edited_codon = ''
        wt_codon = ''
        lst = []
        if len(a_lst) == 1:
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            

        elif len(a_lst) == 2:
            ## edit all 'A'
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            ## edit one 'A'
            wt_codon = ''
            edited_codon = ''
            for idx, nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()     
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])

        else: ## len(a_lst) == 3
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if nt == wt_nt:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()  
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()   
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower() 
                elif idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                           
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[0]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                elif idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                     
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])                    
            wt_codon = ''
            edited_codon = ''
            for idx,nt in enumerate(codon_seq):
                if idx == a_lst[2]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()
                elif idx == a_lst[1]:
                    edited_codon += edited_nt.lower()
                    wt_codon += wt_nt.lower()                          
                else:
                    edited_codon += nt
                    wt_codon += nt
            lst.append([wt_codon,edited_codon])
        
        return lst

    def _edit_all_nt_in_one_codon(self,codon_seq,codon_num):
        a_lst = []
        c_lst = []

        for idx,nt in enumerate(codon_seq):
            if nt == 'T':
                a_lst.append(idx)
            elif nt == 'G':
                c_lst.append(idx)
            else:
                pass
        if len(a_lst) ==0:
            if len(c_lst) == 0:
                return pd.DataFrame()
            else:
                ## C 만 있음
                edited_df = self._edit_every_AorC(c_lst,codon_seq,'G','A')
                edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                edited_df['BaseEditor'] = 'CBE'
        else:
            if len(c_lst) == 0:
                edited_df = self._edit_every_AorC(a_lst,codon_seq,'T','C')
                edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                
                edited_df['BaseEditor'] = 'ABE'
            else: ## 둘 다 있음

                cbe_df = self._edit_every_AorC(c_lst,codon_seq,'G','A')
                cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                cbe_df['BaseEditor'] = 'CBE'

                abe_df = self._edit_every_AorC(a_lst,codon_seq,'T','C')
                abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                abe_df['BaseEditor'] = 'ABE'

                edited_df = pd.concat([cbe_df,abe_df])

        edited_df['Codon_Num'] = codon_num
        return edited_df


    def _edit_partial_nt_in_one_codon(self,codon_seq,codon_num,tup,direction):
        
        final = []
        if direction=='left':
            ## remnant가 왼쪽
            editable_seq_in_codon = codon_seq[tup[0]:]
            remnant =  codon_seq[:tup[0]]
            lst = []
            if len(editable_seq_in_codon) == 1:
                if editable_seq_in_codon == 'T':
                    wt = remnant+'t'
                    edited = remnant+'c'
                    lst.append([wt,edited,'ABE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                elif editable_seq_in_codon == 'G':
                    wt = remnant+'g'
                    edited = remnant+'a'
                    lst.append([wt,edited,'CBE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                else:
                    edited_df = pd.DataFrame()


            elif len(editable_seq_in_codon) == 2:
                a_lst = []
                c_lst = []

                for idx,nt in enumerate(editable_seq_in_codon):
                    if nt == 'T':
                        a_lst.append(idx)
                    elif nt == 'G':
                        c_lst.append(idx)
                    else:
                        pass
                
                if len(a_lst) ==0:
                    if len(c_lst) == 0:
                        edited_df = pd.DataFrame()

                    else:
                        ## C 만 있음
                        edited_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'G','A')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        edited_df['BaseEditor'] = 'CBE'
                else:
                    if len(c_lst) == 0:
                        edited_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'T','C')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        
                        edited_df['BaseEditor'] = 'ABE'
                    else: ## 둘 다 있음
                        cbe_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'G','A')
                        cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                        cbe_df['BaseEditor'] = 'CBE'

                        abe_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'T','C')
                        abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                        abe_df['BaseEditor'] = 'ABE'

                        edited_df = pd.concat([cbe_df,abe_df])
                edited_df.index = [i for i in range(edited_df.shape[0])]
                for idx in edited_df.index:
                    wt = edited_df.loc[idx,'WT']
                    edited = edited_df.loc[idx,'Edited']
                    edited_df.loc[idx,'WT'] = remnant+wt
                    edited_df.loc[idx,'Edited'] = remnant+edited
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)
                        
        else:
            edited_df = pd.DataFrame()
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)            
        
         
            
        if direction=='right':
            ## Remnant r가 오른쪽
            editable_seq_in_codon = codon_seq[:tup[1]-1]
            remnant =  codon_seq[tup[1]-1:]
            lst = []
            if len(editable_seq_in_codon) == 1:
                if editable_seq_in_codon == 'T':
                    wt = 't'+remnant
                    edited = 'c'+remnant
                    lst.append([wt,edited,'ABE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                elif editable_seq_in_codon == 'G':
                    wt = 'g'+remnant
                    edited = 'a'+remnant
                    lst.append([wt,edited,'CBE'])
                    edited_df = pd.DataFrame(lst,columns = ['WT','Edited','BaseEditor'])
                else:
                    edited_df = pd.DataFrame()


            elif len(editable_seq_in_codon) == 2:
                a_lst = []
                c_lst = []

                for idx,nt in enumerate(editable_seq_in_codon):
                    if nt == 'T':
                        a_lst.append(idx)
                    elif nt == 'G':
                        c_lst.append(idx)
                    else:
                        pass
                
                if len(a_lst) ==0:
                    if len(c_lst) == 0:
                        edited_df = pd.DataFrame()

                    else:
                        ## C 만 있음
                        edited_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'G','A')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        edited_df['BaseEditor'] = 'CBE'
                else:
                    if len(c_lst) == 0:
                        edited_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'T','C')
                        edited_df = pd.DataFrame(edited_df,columns = ['WT','Edited'])
                        
                        edited_df['BaseEditor'] = 'ABE'
                    else: ## 둘 다 있음
                        cbe_df = self._edit_every_AorC(c_lst,editable_seq_in_codon,'G','A')
                        cbe_df = pd.DataFrame(cbe_df,columns = ['WT','Edited'])
                        cbe_df['BaseEditor'] = 'CBE'

                        abe_df = self._edit_every_AorC(a_lst,editable_seq_in_codon,'T','C')
                        abe_df = pd.DataFrame(abe_df,columns = ['WT','Edited'])
                        abe_df['BaseEditor'] = 'ABE'

                        edited_df = pd.concat([cbe_df,abe_df])

                edited_df.index = [i for i in range(edited_df.shape[0])]
                for idx in edited_df.index:
                    wt = edited_df.loc[idx,'WT']
                    edited = edited_df.loc[idx,'Edited']
                    edited_df.loc[idx,'WT'] = wt+remnant
                    edited_df.loc[idx,'Edited'] = edited+remnant

            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)
        else:
            edited_df = pd.DataFrame()
            edited_df['Codon_Num'] = codon_num
            final.append(edited_df)   
        final = pd.concat(final)
        return final

    
    def possible_editing_pattern(self,tup,editable_codon,edit_window):
        ## check 'A' or 'C' in edit_target_seq
        num_editable_codon = len(editable_codon)
        lst = []
        if self._check_editable_GorT(edit_window) == 0:
            return pd.DataFrame()
        else:
            pass
        
        for i,idx in enumerate(editable_codon):
            codon_tup = editable_codon[idx]
            codon_seq = codon_tup[1]
            codon_num = codon_tup[0]

            if self._check_editable_GorT(codon_seq) == 0:
                return pd.DataFrame()
            else:
                pass
            if tup[0] > 0:
                if i == 0:
                    eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'left')

                    lst.append(eachdf)
                else:
                    if tup[1] > 0:
                        
                        if i == num_editable_codon-1:
                            eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'right')
                            lst.append(eachdf)
                        else:
                            eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                            lst.append(eachdf)
                    else:
                        eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                        lst.append(eachdf)
                
            else: ## tup[0] <=0
                if tup[1] > 0:
                    if i == num_editable_codon-1:
                        eachdf = self._edit_partial_nt_in_one_codon(codon_seq,codon_num,tup,'right')
                        lst.append(eachdf)
                    else:
                        eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                        lst.append(eachdf)
                else:
                    eachdf = self._edit_all_nt_in_one_codon(codon_seq,codon_num)
                    lst.append(eachdf)

        
        final = pd.concat(lst)
        return final

        

  

        








#%%
def codon_start_number_and_search_area(vardf):
    n = 0
    codon_num = 0
    while n < (vardf.shape[0]):
        remnant = vardf.loc[n,'Remnant_NT']
        trimmed = vardf.loc[n,'trimmed_CDS']

        no_codon = len(trimmed)/3

        if n == 0:
            vardf.loc[n,'CodonStart'] = 1
        else:
            vardf.loc[n,'CodonStart'] = codon_num+1

        if remnant != 0:
            codon_num += no_codon+1
        else:
            codon_num += no_codon

        n+=1
    return vardf

def search_ngg_pam(seq,pamtype):
    ## return N position of NGG
    ## editing window = [-3+NGGidx:NGGidx+5]
    n = 0
    lst = []
    if pamtype == 'NGG':
        pam = 'GG'
    elif pamtype == 'NG':
        pam = 'G'
    else:
        pass
    while n < len(seq)-3+1:
        seq = seq.upper()
        motif = seq[n:n+2]
        if motif == pam:
            lst.append(n-1)
        else:
            pass
        n += 1
    return lst


def search_ccn_pam(seq,pamtype):
    ## return C position of CCN
    ## editing window = [-3+CCNidx:CCNidx+5]
    n = 0
    lst = []
    if pamtype == 'NGG':
        pam = ['CC']
    elif pamtype == 'NG':
        pam = ['AC','TC','GC','CC']
    else:
        pass
    while n < len(seq)-3+1:
        seq = seq.upper()
        motif = seq[n:n+2]
        if motif in pam:
            lst.append(n)
        else:
            pass
        n += 1
    return lst


def codon_indexing_in_sa_ngg(trimmed_cds,start_no,sa_ngg):
    idx = 0
    n = 0
    codon_start = start_no
    dic = {}
    index_correction = sa_ngg.find(trimmed_cds)
    while idx < len(trimmed_cds)-3+1:
        codon = trimmed_cds[idx:idx+3]
        codon_num = codon_start+n

        codon_idx = idx
        dic[codon_idx+index_correction] = (codon_num,codon)

        n += 1
        idx += 3

    return dic


def codon_indexing_in_sa_ccn(trimmed_cds,start_no,sa_ccn):
    idx = 0
    n = 0
    codon_start = start_no
    dic = {}
    index_correction = sa_ccn.find(trimmed_cds)

    while idx < len(trimmed_cds)-3+1:
        codon = trimmed_cds[idx:idx+3]
        codon_num = codon_start+n

        codon_idx = idx
        dic[codon_idx+index_correction] = (codon_num,codon)

        n += 1
        idx += 3

    return dic



def find_editable_codon_index(dic,editable_index):
    first = editable_index[0]
    last  = editable_index[-1]
    lst = []
    for idx in dic:
        if idx >=first-2:
            if idx <=last:
                lst.append(idx)
            else:
                pass
        else:
            pass
    return lst


def editable_index_in_ngg(pamlst):
    dic = {}    
    for pam in pamlst:
        editable_index_lst = [pam-i for i in range(13,18)]
        for idx in editable_index_lst:
            dic[idx] = pam
    return dic



def ngg_patterns(sa_ngg,cds_extended,codon_start,trimmed_cds,pamtype,window_start=3,window_end=7):
    ## window 3 to 7
    pamlst = search_ngg_pam(sa_ngg,pamtype)
    pamlst = [i for i in pamlst if i >=15]  ## NGG index not GG

    if len(pamlst) == 0:
        return pd.DataFrame()
    else:
        pass

    codon_dic = codon_indexing_in_sa_ngg(trimmed_cds,codon_start,sa_ngg)
    saindex = cds_extended.find(sa_ngg)
    
    lst = []
    
    for pam in pamlst:
        true_pam_index = pam+saindex
        gx20 = cds_extended[true_pam_index-20:true_pam_index]
        pam_seq = cds_extended[true_pam_index:true_pam_index+3]

        editable_nt_index = [i for i in range(pam-20+window_start-1,pam-20+window_end)]
        edit_window = gx20[window_start-1:window_end]

        editable_codon_index = find_editable_codon_index(codon_dic,editable_nt_index)

        editable_codon = {}
        for i in editable_codon_index:
            if i in codon_dic:
                editable_codon[i] = codon_dic[i]
            else:
                pass

        edit_target_seq = [editable_codon[i][1] for i in editable_codon ]
        edit_target_seq = "".join(edit_target_seq)

        first_editable_nt_index = editable_nt_index[0]
        first_codon_index = list(editable_codon.keys())[0]

        last_editable_nt_index = editable_nt_index[-1]
        last_codon_index = list(editable_codon.keys())[-1]+2

        pre = first_editable_nt_index-first_codon_index
        post = last_codon_index-last_editable_nt_index

        tup = (pre,post)

        eachdf = EditingPattern_NGG().possible_editing_pattern(tup,editable_codon,edit_window)

        eachdf['GX20'] = gx20
        eachdf['PAM'] = pam_seq

        ## deep input = 24 bp target sequence (1 bp + 20 bp protospacer + PAM + 3 bp) 
        deep_input = cds_extended[true_pam_index-20-1:true_pam_index+3+3]

        eachdf['DeepInput'] = deep_input
        lst.append(eachdf)

    result = pd.concat(lst)

    return result




def ccn_patterns(sa_ccn,cds_extended,codon_start,trimmed_cds,pamtype,window_start=3,window_end=7):
    ## window 3 to 7
    pamlst = search_ccn_pam(sa_ccn,pamtype)
    pamlst = [i for i in pamlst if i <=len(sa_ccn)-20]  ## CCN index "CC"
    
    if len(pamlst) == 0:
        return pd.DataFrame()
    else:
        pass

    codon_dic = codon_indexing_in_sa_ccn(trimmed_cds,codon_start,sa_ccn)
    saindex = cds_extended.find(sa_ccn)    

    lst = []
    
    for pam in pamlst:
        true_pam_index = pam+saindex
        gx20 = cds_extended[true_pam_index+3:true_pam_index+23]
        pam_seq = cds_extended[true_pam_index:true_pam_index+3]

        editable_nt_index = [i for i in range(pam+20-window_end+3,pam+20-window_start+3+1)]
        edit_window = gx20[20-window_end:20-window_start+1]

        editable_codon_index = find_editable_codon_index(codon_dic,editable_nt_index)

        editable_codon = {}
        for i in editable_codon_index:
            if i in codon_dic:
                editable_codon[i] = codon_dic[i]
            else:
                pass

        edit_target_seq = [editable_codon[i][1] for i in editable_codon ]
        edit_target_seq = "".join(edit_target_seq)

        first_editable_nt_index = editable_nt_index[0]
        first_codon_index = list(editable_codon.keys())[0]

        last_editable_nt_index = editable_nt_index[-1]
        last_codon_index = list(editable_codon.keys())[-1]+2

        pre = first_editable_nt_index-first_codon_index
        post = last_codon_index-last_editable_nt_index

        tup = (pre,post)

        eachdf = EditingPattern_CCN().possible_editing_pattern(tup,editable_codon,edit_window)
    
        eachdf['GX20'] = gx20
        eachdf['PAM'] = pam_seq

        ## deep input = 24 bp target sequence (1 bp + 20 bp protospacer + PAM + 3 bp) 
        deep_input = cds_extended[true_pam_index-20-1:true_pam_index+3+3]

        eachdf['DeepInput'] = deep_input
        lst.append(eachdf)

    result = pd.concat(lst)

    return result




def make_editable_pattern_in_each_cds(sa_ngg,sa_ccn,cds_extended,codon_start,trimmed_cds,pamtype,window_start,window_end):

    lst = []
    ngg = ngg_patterns(sa_ngg,cds_extended,codon_start,trimmed_cds,pamtype,window_start,window_end)
    ccn = ccn_patterns(sa_ccn,cds_extended,codon_start,trimmed_cds,pamtype,window_start,window_end)
    
    pam = pamtype
    rt_pam = General().rt_sequence(pamtype)

    ngg['PAM'] = pam
    ccn['PAM'] = rt_pam

    lst.append(ngg)
    lst.append(ccn)

    if len(lst) == 0:
        return pd.DataFrame()
    else:
        final = pd.concat(lst)
        return final

    
    
    


def make_all_possible_editing_result(eachvar,var_num,pamtype,window_start,window_end):
    ##ccn add
    lst = []
    for idx in eachvar.index:
        
        trimmed_cds = eachvar.loc[idx,'trimmed_CDS']
        cds_extended    = eachvar.loc[idx,'CDS_extended']

        trimmed_start = cds_extended.find(trimmed_cds)

        sa_ngg = cds_extended[trimmed_start:trimmed_start+len(trimmed_cds)+20]
        sa_ccn = cds_extended[trimmed_start-20:trimmed_start+len(trimmed_cds)]

        
        codon_start = eachvar.loc[idx,'CodonStart']

        result = make_editable_pattern_in_each_cds(sa_ngg,sa_ccn,cds_extended,codon_start,trimmed_cds,pamtype,window_start,window_end)

        lst.append(result)


    final = pd.concat(lst)
    final['Variants'] = var_num
    
    return final


#%%

def pe_editable_pattern(eachdf,pamtype,window_start,window_end):
    t1 = datetime.now()    
    dic = extract_cds_from_all_variants().cds_extractor(eachdf)

    lst = []
    for var_num in dic:
        eachvar = dic[var_num]
        if eachvar.shape[0] == 0:
            return pd.DataFrame()

        eachvar = codon_start_number_and_search_area(eachvar)
        eachvar = make_all_possible_editing_result(eachvar,var_num,pamtype,window_start,window_end)
        if eachvar.shape[0] == 0:
            pass
        else:
            lst.append(eachvar)

    t2 = datetime.now()
    if len(lst) == 0:
        final = pd.DataFrame()
    else:
        final = pd.concat(lst)
    
    print(t2-t1)
    return final

#%%
def main(genelst,pamtype,window_start,window_end,output_name):
    df = pd.read_csv('EssentialData/CCDS.current.txt',sep='\t')
    df = df[df['cds_locations']!='-']
    abelst = []
    cbelst = []



    with ProcessPoolExecutor(max_workers=None) as executor:
        futs = []
        for gene in genelst:
            if gene not in list(df['gene']):
                pass
            else:
                eachdf = df[df['gene']==gene]
                fut = executor.submit(pe_editable_pattern,eachdf,pamtype,window_start,window_end)
                futs.append((gene,fut))

        
        for tup in futs:
            gene = tup[0]
            result = tup[1].result()
            if result.shape[0] == 0:
                pass
            else:
                result['gene'] = gene

                cbe = result[result['BaseEditor']=='CBE']
                abe = result[result['BaseEditor']=='ABE']

                cbelst.append(cbe)
                abelst.append(abe)
    
    abe = pd.concat(abelst)
    cbe = pd.concat(cbelst)

    abe.to_csv('Output/{}_ABE_Tiling_Library.txt'.format(output_name),sep='\t',index=None)
    cbe.to_csv('Output/{}_CBE_Tiling_Library.txt'.format(output_name),sep='\t',index=None)

    return
# %%
if __name__ == '__main__':
    pamtype = 'NGG'
    window_start = 3
    window_end = 7

    output_name = input('ouput_name : ')


    
    with open('Input/input.txt','r') as inp:
        genelst = inp.readlines()
        genelst = [i.rstrip() for i in genelst if i.rstrip() !='gene']
        
    
    main(genelst,pamtype,window_start,window_end,output_name)

#%%

# %%
