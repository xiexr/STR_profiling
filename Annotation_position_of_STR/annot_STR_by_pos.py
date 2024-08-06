import re
import sys
from rice_IRGSP_gff import gffInfo
import collections
 
region_seq = ["CDS", "five_prime_UTR", "three_prime_UTR", "exon", "pos"]


class gffJudge(gffInfo):
    def __init__(self, var_fn ):#gff_info_dic
        self.var_fn = var_fn 
        #self.gff_info_dic = gff_info_dic 
        # get the gff_info dict
        super(gffJudge, self).__init__()
        # judge
        self.judgePos()

    def judgePos(self):
        pos_dict = collections.OrderedDict() 
        var_dict = collections.OrderedDict()
        for l in open(self.var_fn, "r"):
            var = l.rstrip()
            onlynum = re.search("\d+", var).group(0)
            chrom = "chr%s" % onlynum[:2]  # chr01
            pos=str(int(onlynum[2:]))
            pos_dict.setdefault(chrom, []).append(int(pos))#{chr12:[1266,1326,1383,1396,1512]}
            var_dict.setdefault(chrom, []).append(var)#{chr12:[vg0100001266,vg0100002369]}
            #print chrom, pos

        # annot the pos of each chrom
        for chrom, pos_lst in pos_dict.items():
            var_lst = var_dict[chrom] #[vg0100002722,vg010023132156,vg23132456]
            self.judgePosEach(chrom, var_lst, pos_lst)
            # hit_region_dict = self.judgePosEach(chrom, var_lst, pos_lst) #chr01,[vg0100002722,vg123456456,vg54596463],[2722,546465,123456456,45485]
    def judgePosEach(self, chrom, varQuery_lst, posQuery_lst):
        
        anote_dic = {}   # a container to temporary record the annotation info
        for i in posQuery_lst:
            anote_dic.setdefault(i, {}) #{2722:{},546465:{}}
        for j in range(len(posQuery_lst)):#遍历所有STR位置
            pos = posQuery_lst[j] #2722
            var = varQuery_lst[j] #vg0100002722
            for g in range(len(self.gene_order_dic[chrom])):#遍历所有基因
                gene_id = self.gene_order_dic[chrom][g]
                gene_symb = self.gene_symb_dic[gene_id][1]
                gene_descript = self.gene_symb_dic[gene_id][0]
                if len(gene_symb) == 0:
                    gene_symb = gene_id
                gene_pos = self.gene_pos_dic[gene_id]# ['chr12', 5685138, 5689991, '-']
                bef_gene = self.bef_gene_dic[chrom][gene_id]#相邻基因,该基因前一个
                # trans_dic = self.gff_dic[chrom][gene_id]
                # alltranspos = [(i, int(trans_dic.get(i).get('pos')[0][1]) - int(trans_dic.get(i).get('pos')[0][0])) for
                #                i in trans_dic]
                # sortalltranspos = sorted(alltranspos, key=lambda x: x[1], reverse=True)
                # longesttransposid = sortalltranspos[0][0]
                # info_dic = trans_dic[longesttransposid]
                # trans_pos = info_dic["pos"][0]
                # bef_trans_dic=self.gff_dic[chrom][bef_gene]
                # bef_alltranspos=[(i, int(bef_trans_dic.get(i).get('pos')[0][1]) - int(bef_trans_dic.get(i).get('pos')[0][0])) for
                #                i in bef_trans_dic]
                # bef_sortalltranspos = sorted(bef_alltranspos, key=lambda x: x[1], reverse=True)
                # bef_longesttransposid = bef_sortalltranspos[0][0]
                # bef_info_dic=bef_trans_dic[bef_longesttransposid]
                # bef_trans_pos = bef_info_dic["pos"][0]
                if pos > int(gene_pos[2]):#大于该基因的末尾
                    if g == len(self.gene_order_dic[chrom])-1:
                        region = "intergenic"
                        loc_region = "%s__" % gene_id
                        symb_region = "%s__" % gene_symb
                        descript = "%s__:" % gene_descript
                        trans_dic = self.gff_dic[chrom][gene_id]
                        alltranspos = [(i, int(trans_dic.get(i).get('pos')[0][1]) - int(trans_dic.get(i).get('pos')[0][0])) for i in trans_dic]
                        sortalltranspos = sorted(alltranspos, key=lambda x: x[1], reverse=True)
                        longesttransposid = sortalltranspos[0][0]
                        info_dic = trans_dic[longesttransposid]
                        trans_pos = info_dic["pos"][0]
                        right_gene_pos=int(trans_pos[-1])
                        if gene_pos[3] == "-":
                            left_stream = "upstream"
                            left_distance = pos - right_gene_pos
                        else:
                            left_stream = "downstream"
                            left_distance = pos - right_gene_pos
                        right_gene = "" 
                        right_stream = ""
                        right_distance = ""
                        print( "%s\t%s\t%s\t%s\t%s\t%s:%s_%s;%s:%s_%s\t%s" % (var, chrom, pos, loc_region, symb_region, gene_id, left_stream, left_distance, right_gene, right_stream, right_distance, descript))
                    #% s: % s_ % s; % s: % s_ % s,gene_id,left_stream, left_distance, right_gene, right_stream, right_distance
                    continue
                else:
                    if pos >= int(gene_pos[1]):#大于开头小于末尾
                        gene_strand = gene_pos[3]    # the direct of transcript == gene
                        trans_dic = self.gff_dic[chrom][gene_id] 
                        annot_lst = []
                        alltranspos=[(i,int(trans_dic.get(i).get('pos')[0][1])-int(trans_dic.get(i).get('pos')[0][0])) for i in trans_dic]
                        sortalltranspos=sorted(alltranspos,key=lambda x:x[1],reverse=True)
                        longesttransposid=sortalltranspos[0][0]
                        info_dic = trans_dic[longesttransposid]
                        # get the trans pos
                        trans_pos = info_dic["pos"][0]
                        if gene_strand=='+':
                            if pos < int(trans_pos[0]):
                                region = "%s:up_%s" % (longesttransposid, int(trans_pos[0]) - pos)
                            elif pos > int(trans_pos[1]):
                                region = "%s:down_%s" % (longesttransposid, pos - int(trans_pos[1]))
                            else:
                                hit_region_lst = []
                                for t in info_dic.keys():
                                    pos_lst = info_dic[t]
                                    for each_lst in pos_lst:
                                        start = int(each_lst[0])
                                        end = int(each_lst[1])
                                        if start<=pos<=end:
                                            hit_region_lst.append(t)
                                for k in region_seq:
                                    if k in hit_region_lst:
                                        region = k
                                        break
                                if region == "pos":
                                    region = "%s:intron:%s" % (longesttransposid,pos-int(trans_pos[0]))
                                else:
                                    region = "%s:%s:%s" % (longesttransposid, region,pos-int(trans_pos[0]))
                        else:
                            if pos < int(trans_pos[0]):
                                region = "%s:down_%s" % (longesttransposid, int(trans_pos[0]) - pos)
                            elif pos > int(trans_pos[1]):
                                region = "%s:up_%s" % (longesttransposid, pos - int(trans_pos[1]))
                            else:
                                hit_region_lst = []
                                for t in info_dic.keys():
                                    pos_lst = info_dic[t]
                                    for each_lst in pos_lst:
                                        start = int(each_lst[0])
                                        end = int(each_lst[1])
                                        if start <= pos <= end:
                                            hit_region_lst.append(t)
                                for k in region_seq:
                                    if k in hit_region_lst:
                                        region = k
                                        break
                                if region == "pos":
                                    region = "%s:intron:%s" % (longesttransposid,int(trans_pos[1])-pos)
                                else:
                                    region = "%s:%s:%s" % (longesttransposid, region,int(trans_pos[1])-pos)
                        annot_lst.append(region)

                        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (var, chrom, pos, gene_id, gene_symb, ";".join(annot_lst), gene_descript))
                    else:
                        region = "intergenic"
                        loc_region = "%s__%s" % (bef_gene, gene_id)
                        if len(bef_gene) == 0:
                            bef_gene_symb = ""
                            bef_gene_descript = ""
                        else:
                            bef_gene_symb = self.gene_symb_dic[bef_gene][1]
                            bef_gene_descript = self.gene_symb_dic[bef_gene][0]
                            if len(bef_gene_symb) == 0:
                                bef_gene_symb = bef_gene
                        symb_region = "%s__%s" % (bef_gene_symb, gene_symb)
                        descript = "%s__%s" % (bef_gene_descript, gene_descript)

                        # the distance to the left and right gene
                        if len(bef_gene) != 0:
                            bef_trans_dic = self.gff_dic[chrom][bef_gene]
                            bef_alltranspos = [
                                (i, int(bef_trans_dic.get(i).get('pos')[0][1]) - int(bef_trans_dic.get(i).get('pos')[0][0])) for
                                i in bef_trans_dic]
                            bef_sortalltranspos = sorted(bef_alltranspos, key=lambda x: x[1], reverse=True)
                            bef_longesttransposid = bef_sortalltranspos[0][0]
                            bef_info_dic = bef_trans_dic[bef_longesttransposid]
                            bef_trans_pos = bef_info_dic["pos"][0]
                            left_gene_pos = int(bef_trans_pos[-1])
                            bef_gene_pos=self.gene_pos_dic[bef_gene]
                            if bef_gene_pos[3] == "-":
                                left_stream = "upstream"
                                left_distance = pos - left_gene_pos
                            else:
                                left_stream = "downstream"
                                left_distance = pos - left_gene_pos
                        else:
                            left_stream = ""
                            left_distance = ""
                        trans_dic = self.gff_dic[chrom][gene_id]
                        alltranspos = [
                            (i, int(trans_dic.get(i).get('pos')[0][1]) - int(trans_dic.get(i).get('pos')[0][0])) for
                            i in trans_dic]
                        sortalltranspos = sorted(alltranspos, key=lambda x: x[1], reverse=True)
                        longesttransposid = sortalltranspos[0][0]
                        info_dic = trans_dic[longesttransposid]
                        trans_pos = info_dic["pos"][0]
                        right_gene_pos = int(trans_pos[0])
                        if gene_pos[3] == "+":
                            right_stream = "upstream"
                            right_distance = right_gene_pos - pos
                        else:
                            right_stream = "downstream"
                            right_distance = right_gene_pos - pos
                        print ('%s\t%s\t%s\t%s\t%s\t%s:%s_%s;%s:%s_%s\t%s' % (var, chrom, pos, loc_region, symb_region, bef_gene, left_stream, left_distance, gene_id, right_stream, right_distance, descript))
                    break 


################################################################################################
(xx, var_fn)=sys.argv

res = gffJudge(var_fn)
 #######