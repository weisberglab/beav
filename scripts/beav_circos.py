#!/usr/bin/python3
import argparse
import sys, os
sys.path.insert(0, f'{os.environ["BEAV_DIR"]}/software/PyCirclize/src/')


from pycirclize import Circos
from pycirclize.parser import Genbank
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import beav_oncogenes

def get_base_file_name(file_name):
    """
    INPUT: Genbank file name or path to the file
    OUTPUT: Base name 
    """
    if '.' in file_name:
        file_name = file_name.split('.')[0]
    
    if '/' in file_name:
        file_name = file_name.split('/')[-1]
    
    return file_name

def all_contig_circos(gbk_file):
    """
    INPUT: Genbank file containing beav annotation
    OUTPUT: Circos plot containing feature distribution in differnt contigs
    """
    # Load GBK file
    gbk = Genbank(gbk_file)

    # Get contig genome seqid & size, features dict
    seqid2size = gbk.get_seqid2size()
    seqid2features = gbk.get_seqid2features(feature_type=None)


    circos = Circos(seqid2size, space=1, start=15, end=345)
    #main_kwargs = dict(ha="right", va="top")
    #circos.text(f"{get_base_file_name(gbk_file)}\n{len(circos.sectors)} contig(s)", r=12, size=8, **main_kwargs)

    # Loop through the contigs
    contig_i = 0
    for sector in circos.sectors:
        # plot sector labels
        kwargs = dict(color="navy")
        sector.text(f"{sector.name}", orientation='vertical', r=110, size=10, **kwargs)
        # Plot outer track
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fc="lightgrey")
        major_interval = 1000000
        minor_interval = int(major_interval / 100)
        if sector.size > minor_interval:
            outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000000:.0f} Mb")
            outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
    
        f_cds_track = sector.add_track((90, 95), r_pad_ratio=0.1)
        r_cds_track = sector.add_track((85, 90), r_pad_ratio=0.1)
        beav_track = sector.add_track((80, 85), r_pad_ratio=0.1)
        extra_feature =  sector.add_track((72, 77))
        gc_content_track = sector.add_track((67, 72))
        gc_skew_track = sector.add_track((62, 67))
        
        # Plot forward/reverse CDS, ICE, origin, beav tracks
        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                if feature.strand == 1:
                    f_cds_track.genomic_features([feature], fc="tomato")
                else:
                    r_cds_track.genomic_features([feature], fc="skyblue")
            elif feature.type == "mobile_element":
                    if 'mobile_element_type' in feature.qualifiers.keys():
                        if feature.qualifiers['mobile_element_type'][0] == 'phage':
                            beav_track.genomic_features([feature], fc='royalblue') # Pahge
                        elif feature.qualifiers['mobile_element_type'][0] == 'integrative element':
                            beav_track.genomic_features([feature], fc="turquoise") # Non phage mobile element
            elif feature.type == "rRNA":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=True, label_size=6, labels=[''], label_orientation="vertical", line_kws={'ec':'yellowgreen'}) #rRNA
            elif feature.type == "tRNA":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=True, label_size=6, labels=[''], label_orientation="vertical", line_kws={'ec':'orange'}) #tRNA

                
            if feature.type == "CDS":
                
                if "note" in feature.qualifiers.keys(): 
                    if 'MacSy' in feature.qualifiers['note'][0]:
                        sec_sys = feature.qualifiers['note'][0].split(': ')[1]
                        beav_track.genomic_features([feature], fc="darkorange") # Secretion Systems
                        if 'T4SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # T4SS
                        elif 'T3SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'plum'}, text_kws={'color':'darkred'}) # T3SS
                        elif 'T6SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'slateblue'}, text_kws={'color':'darkred'}) # T6SS
                    elif 'SMASH' in feature.qualifiers['note'][0]:
                        beav_track.genomic_features([feature], fc="forestgreen") #Secondary metabolite
            
            if feature.type == "rep_origin":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['origin'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # Origin of replication
            if feature.type == "CDS":
                if 'gene' in feature.qualifiers.keys():
                    if feature.qualifiers['gene'][0] in ['repA', 'repB', 'repC']:
                        extra_feature.genomic_features([feature], fc='black',  plotstyle="arrow")
                                            

        # Plot GC skew
        pos_list, gc_skews = gbk.calc_gc_skew(seq=str(gbk.records[contig_i].seq))
        positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
        negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
        abs_max_gc_skew = np.max(np.abs(gc_skews))
        vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
        gc_skew_track.fill_between(
            pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
        )
        gc_skew_track.fill_between(
            pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
        )


        # update contig index
        contig_i += 1


    text_common_kws = dict(ha="center", va="center", size=6)
    circos.text("CDS +strand", r=93, color="dimgrey", **text_common_kws)
    circos.text("CDS -strand", r=87, color="dimgrey", **text_common_kws)
    circos.text("Beav", r=82, color="dimgrey", **text_common_kws)
    circos.text("repABC", r=74, color="dimgrey", **text_common_kws)
    circos.text("RNA", r=68, color="dimgrey", **text_common_kws)
    circos.text("GC skew", r=63, color="dimgrey", **text_common_kws)

    fig = circos.plotfig()
    _ = circos.ax.legend(
        handles=[
            Patch(color="darkorange", label="Secretion Systems"),
            Line2D([], [], color="darkred", label="T4SS", marker='_', ls='None'),
            Line2D([], [], color="plum", label="T3SS", marker='_', ls='None'),
            Line2D([], [], color="slateblue", label="T6SS", marker='_', ls='None'),
            Patch(color="royalblue", label="Phages"),
            Patch(color="turquoise", label="ICEs"),
            Patch(color="forestgreen", label="Secondary Metabolites"),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
            Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
            Line2D([], [], color="yellowgreen", label="rRNA", marker="_", ms=6, ls="None"),
            Line2D([], [], color="orange", label="tRNA", marker="_", ms=6, ls="None")
        ],
        bbox_to_anchor=(0.5, 0.45),
        loc="center",
        #ncols=2,
        fontsize=6
    )

    fig.savefig(f'{get_base_file_name(gbk_file)}.circo.png')
    print(f'Image saved: {get_base_file_name(gbk_file)}.circo.png')

    pass


def splice_genbank_contig(gbk_file, contig_id):
    """
    INPUT:
        gbk_file: Genbank file containing beav annotation and has multiple contigs
        contig_id: List of contigs
    OUTPUT: A spliced genbank file that has only the contig of interest
    """
    with open(gbk_file, 'r') as f:
        gbk_raw = f.read()
    
    gbk_contigs = gbk_raw.strip().split('//\n')
    
    contig_of_interest = []
    for i in gbk_contigs:
        for j in contig_id:
            if j in i.split('\n', 1)[0].split():
                contig_of_interest.append(i)

    with open(f'{get_base_file_name(gbk_file)}.oncogenic.gbk', 'w') as f:
        for i in range(len(contig_of_interest)):
            f.writelines(contig_of_interest[i])
            if i < len(contig_of_interest):
                f.writelines('//\n')
            
    
    print(f'{get_base_file_name(gbk_file)}.oncogenic.gbk')

    return f'{get_base_file_name(gbk_file)}.oncogenic.gbk'


def oncogenic_circos(gbk_file):
    """
    INPUT: 
    gbk_file: Single contig genbank file containing beav annotation
    OUTPUT: Circos plot containing oncogenic plasmid feature distribution of that contig
    """

    #Load GBK file
    gbk = Genbank(gbk_file)

    #Get features
    seqid2size = gbk.get_seqid2size()
    seqid2features = gbk.get_seqid2features(feature_type=None)
    
    #Define sector for single contig
    circos = Circos(seqid2size, space=3)
    #main_kwargs = dict(ha="center", va="top") 
    #circos.text(f"{get_base_file_name(gbk_file)}\nTi/Ri plasmid", r=5, size=8, **main_kwargs)

    # Extract CDS gene labels
    # List containing CDS product labels if gene names are not present
    product_of_interest = []

    for sector in circos.sectors:
        # plot sector labels
        kwargs = dict(color="navy")
        sector.text(f"{sector.name}", orientation='vertical', r=110, size=10, **kwargs)

        # add tracks
        cds_track = sector.add_track((95, 100))
        scale_track = sector.add_track((95, 95))
        cds_track.axis(fc="#EEEEEE", ec="none") 

        # Get sector specific labels
        ##print(sector, seqid2features[sector.name])

        sector_pos_list, sector_labels = [], []
        for feat in seqid2features[sector.name]:
            start, end = int(str(feat.location.end)), int(str(feat.location.start))
            pos = (start + end) / 2
            label = feat.qualifiers.get("gene", [""])[0]
            product = feat.qualifiers.get("product", [""])[0]
            if label == "":
                if product in product_of_interest:
                    label = product
            if product == "" or product.startswith("hypothetical"):
                continue
            if len(label) > 20:
                label = label[:20] + "..."
            sector_pos_list.append(pos)
            sector_labels.append(label)

        # plot oncogene specific features
        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                cds_track.genomic_features([feature], fc="lightgrey", plotstyle="arrow")

            if 'gene' in feature.qualifiers.keys():
                if feature.qualifiers['gene'][0] in beav_oncogenes.vir_dict.keys():
                    cds_track.genomic_features([feature], fc="olive", plotstyle="arrow") #T-DNA transfer
                elif feature.qualifiers['gene'][0] in beav_oncogenes.tra_dict.keys():
                    cds_track.genomic_features([feature], fc="indigo", plotstyle="arrow") #tra genes
                elif feature.qualifiers['gene'][0] in beav_oncogenes.trb_dict.keys():
                    cds_track.genomic_features([feature], fc="purple", plotstyle="arrow") #trb genes
                elif feature.qualifiers['gene'][0] in beav_oncogenes.oncogene_dict.keys():
                    cds_track.genomic_features([feature], fc="orange", plotstyle="arrow") #T-DNA/Oncogene            
                elif feature.qualifiers['gene'][0] in beav_oncogenes.rep_dict.keys():
                    cds_track.genomic_features([feature], fc="salmon", plotstyle="arrow") #repABC
                elif feature.qualifiers['gene'][0] in beav_oncogenes.opine_synth_dict.keys():
                    cds_track.genomic_features([feature], fc="darkgreen", plotstyle="arrow") #Opine synthase genes
                elif feature.qualifiers['gene'][0] in beav_oncogenes.opine_tracat_dict.keys():
                    cds_track.genomic_features([feature], fc="turquoise", plotstyle="arrow") #Opine transport/catabolism genes
                elif feature.qualifiers['gene'][0] in beav_oncogenes.agrocinopine_dict.keys():
                    cds_track.genomic_features([feature], fc="steelblue", plotstyle="arrow") #Agrocinopine transport/catabolism genes
                
            if 'product' in feature.qualifiers.keys():
                if feature.qualifiers['product'][0] in beav_oncogenes.oncogene_list:
                    product_of_interest.append(feature.qualifiers['product'][0]) 
                    cds_track.genomic_features([feature], fc="orange", plotstyle="arrow") #T-DNA/Oncogene
                elif feature.qualifiers['product'][0] in beav_oncogenes.rep_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="salmon", plotstyle="arrow") #repABC
                elif feature.qualifiers['product'][0] in beav_oncogenes.opine_synth_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="darkgreen", plotstyle="arrow") #Opine synthase genes
                elif feature.qualifiers['product'][0] in beav_oncogenes.opine_tracat_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="turquoise", plotstyle="arrow") #Opine transport/catabolism genes
                elif feature.qualifiers['product'][0] in beav_oncogenes.agrocinopine_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="steelblue", plotstyle="arrow") #Agrocinopine transport/catabolism genes
                elif 'transposase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Transposase
                elif 'integrase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Integrase
                elif 'recombinase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Recombinase
                elif 'ABC transporter' in feature.qualifiers['product'][0]:
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'darkgreen'}) #ABC Transporter
            
            if 'note' in feature.qualifiers.keys():
                if 'origin_of_replication' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['origin'], label_orientation="vertical") # Origin of replication
                elif 'virbox' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['virbox'], label_orientation="vertical", line_kws={'ec':'navy'}, text_kws={'color':'navy'}) # virbox
                elif 'trabox' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['trabox'], label_orientation="vertical", line_kws={'ec':'navy'}, text_kws={'color':'navy'}) # trabox
                elif 'T-DNA_right_border' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['T-DNA Right'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # t-dna right border
                elif 'T-DNA_left_border' in feature.qualifiers['note']:
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['T-DNA Left'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # t-dna left border
                elif 'integrase' in ' '.join(feature.qualifiers['note']).lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Integrase

        # Plot CDS product labels on outer position
        cds_track.xticks(
            sector_pos_list,
            sector_labels,
            label_orientation="vertical",
            show_bottom_line=True,
            outer=True,
            label_size=6,
            line_kws=dict(ec="grey"),
        )
        # Plot xticks & intervals on inner position
        scale_track.xticks_by_interval(
            interval=25000,
            label_size=6,
            outer=False,
            show_bottom_line=True,
            label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
            label_orientation="vertical",
            line_kws=dict(ec="grey"),
        )


    fig = circos.plotfig()
    _ = circos.ax.legend(
        handles=[
            Line2D([], [], color="lightgrey", label="CDS +strand", marker=">", ms=6, ls="None"),
            Line2D([], [], color="lightgrey", label="CDS -strand", marker="<", ms=6, ls="None"),
            Patch(color="olive", label="T-DNA transfer"),
            Patch(color="orange", label="T-DNA/Oncogenes"),
            Line2D([], [], color="darkgreen", label="ABC Transporter", marker="_", ms=6, ls="None"),
            Patch(color="indigo", label="tra"),
            Patch(color="purple", label="trb"),
            Patch(color="salmon", label="repABC"),
            Patch(color="darkgreen", label="Opine synthase "),
            Patch(color="turquoise", label="Opine transport"),
            Patch(color="steelblue", label="Agrocinopine transport"),
            Patch(color="gray", label="Transposase/Recombinase"), 
        ],
        bbox_to_anchor=(0.5, 0.45),
        loc="center",
        #ncols=2,
        fontsize=6
    )
    fig.savefig(f'{get_base_file_name(gbk_file)}.oncogenes.png')
    print(f'Image saved: {get_base_file_name(gbk_file)}.oncogenes.png')

    pass

if __name__ == "__main__":
    ## Define arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', type=str, required=True)
    parser.add_argument('--contigs', '-c', type=str, nargs='*', required=False)

    # Parse the arguments
    args = parser.parse_args()

    # Generate only genome circos 
    # Command `python3 beav_circos.py --input input.gbk`
    if not args.contigs:
        print('All contig circos only ...')
        all_contig_circos(args.input)

    # Generate genome and oncogene circos
    # Command: `python3 beav_circos.py --input input.gbk --contig contig_of_interest`
    if args.contigs:
        # Generate genome circos
        all_contig_circos(args.input) 
        # Subset contig_of_interest and create contig_id.gbk file
        contig_id_gbk = splice_genbank_contig(args.input, args.contigs)
        # Generate circos for contig_of_interest
        oncogenic_circos(contig_id_gbk)
