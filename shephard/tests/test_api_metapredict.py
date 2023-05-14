import shephard
from shephard.apis import metapredict as meta_api
from shephard.apis import uniprot
import pytest
from shephard.exceptions import ProteinException

import numpy as np


IDR_seqs = {'O00401-1-14': 'MSSVQQQPPPPRRV',
            'O00401-141-217': 'RQRKSEKRRDPPNGPNLPMATVDIKNPEITTNRFYGPQVNNISHTKEKKKGKAKKKRLTKADIGTPSNFQHIGHVGW',
            'O00401-270-505': 'NELRRQAPPPPPPSRGGPPPPPPPPHNSGPPPPPARGRGAPPPPPSRAPTAAPPPPPPSRPSVAVPPPPPNRMYPPPPPALPSSAPSGPPPPPPSVLGVGPVAPPPPPPPPPPPGPPPPPGLPSDGDHQVPTTAGNKAALLDQIREGAQLKKVEQNSRPVSCSGRDALLDQIRQGIQLKSVADGQESTPPTPAPTSGIVGALMEVMQKRSKAIHSSDEDEDEDDEEDFEDDDEWED',
            'O00470-1-113': 'MAQRYDDLPHYGGMDGVGIPSTMYGDPHAARSMQPVHHLNHGPPLHSHQYPHTAHTNAMAPSMGSSVNDALKRDKDAIYGHPLFPLLALIFEKCELATCTPREPGVAGGDVCS',
            'O00470-184-279': 'DLVIDDREGGSKSDSEDITRSANLTDQPSWNRDHDDTASTRSGGTPGPSSGGHTSHSGDNSSEQGDGLDNSVASPSTGDDDDPDKDKKRHKKRGIF',
            'O00470-333-390': 'PMIDQSNRAVSQGTPYNPDGQPMGGFVMDGQQHMGIRAPGPMSGMGMNMGMEGQWHYM',
            'O00472-139-204': 'TQAEEESRNRSTKVIKPGGPYVGKRVQIRKAPQAVSDTVPERKRSTPMNPANTIRKTHSSSTISQR',
            'O00472-285-521': 'SVLSRKLNPSQNAAGTSRSESPVCSSRDAVSSPQKRLLDSEFIDPLMNKKARISHLTNRVPPTLNGHLNPTSEKSAAGLPLPPAAAAIPTPPPLPSTYLPISHPPQIVNSNSNSPSTPEGRGTQDLPVDSFSQNDSIYEDQQDKYTSRTSLETLPPGSVLLKCPKPMEENHSMSHKKSKKKSKKHKEKDQIKKHDIETIEEKEEDLKREEEIAKLNNSSPNSSGGVKEDCTASMEPS',
            'O00499-273-525': 'QHGSNTFTVKAQPSDNAPAKGNKSPSPPDGSPAATPEIRVNHEPEPAGGATPGATLPKSPSQLRKGPPVPPPPKHTPSKEVKQEQILSLFEDTFVPEISVTTPSQFEAPGPFSEQASLLDLDFDPLPPVTSPVKAPTPSGQSIPWDLWEPTESPAGSLPSGEPSAAEGTFAVSWPSQTAEPGPAQPAEASEVAGGTQPAAGAQEPGETAASEAASSSLPAVVVETFPATVNGTVEGGSGAGRLDLPPGFMFKV',
            'O00629-1-61': 'MADNEKLDNQRLKNFKNKGRDLETMRRQRNEVVVELRKNKRDEHLLKRRNVPHEDICEDSD',
            'O00629-489-521': 'DIDEDPSLVPEAIQGGTFGFNSSANVPTEGFQF',
            'O00712-186-420': 'EQDSGQSGSPSHNDPAKNPPGYLEDSFVKSGVFNVSELVRVSRTPITQGTGVNFPIGEIPSQPYYHDMNSGVNLQRSLSSPPSSKRPKTISIDENMEPSPTGDFYPSPSSPAAGSRTWHERDQDMSSPTTMKKPEKPLFSSASPQDSSPRLSTFPQHHHPGIPGVAHSVISTRTPPPPSPLPFPTQAILPPAPSSYFSHPTIRYPPHLNPQDTLKNYVPSYDPSSPQTSQSWYLG',
            'O00716-1-178': 'MRKGIQPALEQYLVTAGGGEGAAVVAAAAAASMDKRALLASPGFAAAAAAAAAPGAYIQILTTNTSTTSCSSSLQSGAVAAGPLLPSAPGAEQTAGSLLYTTPHGPSSRAGLLQQPPALGRGGSGGGGGPPAKRRLELGESGHQYLSDGLKTPKGKGRAALRSPDSPKTPKSPSEKTR',
            'O00716-351-465': 'CPEETETHSPMKTNNQDHNGNIPKPASKDLASTNSGHSDCSVSMGNLSPLASPANLLQQTEDQIPSNLEGPFVNLLPPLLQEDYLLSLGEEEGISDLFDAYDLEKLPLVEDFMCS',
            'O14786-1-24': 'MERGLPLLCAVLALVLAPAGAFRN',
            'O14786-809-855': 'CAKPADLDKKNPEIKIDETGSTPGYEGEGEGDKNISRKPGNVLKTLD',
            'O14786-887-923': 'GMSERNLSALENYNFELVDGVKLKKDKLNTQSTYSEA',
            'Q9UJX3-111-132': 'SKTSKVRPSTGNSASTPQSQCL',
            'Q9UJX3-540-599': 'PNDQKSLEGMQKMEKEESPTDATQEEDVDDMEGSGEEGDLEGSDSEAAQWADQEQWFGMQ'}

folded_seqs = {'O00401-15-140': 'TNVGSLLLTPQENESLFTFLGKKCVTMSSAVVQLYAADRNCMWSKKCSGVACLVKDNPQRSYFLRIFDIKDGKLLWEQELYNNFVYNSPRGYFHTFAGDTCQVALNFANEEEAKKFRKAVTDLLGR',
               'O00401-218-269': 'DPNTGFDLNNLDPELKNLFDMCGISEAQLKDRETSKVIYDFIEKTGGVEAVK',
               'O00470-114-183': 'SESFNEDIAVFAKQIRAEKPLFSSNPELDNLMIQAIQVLRFHLLELEKVHELCDNFCHRYISCLKGKMPI',
               'O00470-280-332': 'PKVATNIMRAWLFQHLTHPYPSEEQKKQLAQDTGLTILQVNNWFINARRRIVQ',
               'O00472-1-138': 'MAAGGTGGLREEQRYGLSCGRLGQDNITVLHVKLTETAIRALETYQSHKNLIPFRPSIQFQGLHGLVKIPKNDPLNEVHNFNFYLSNVGKDNPQGSFDCIQQTFSSSGASQLNCLGFIQDKITVCATNDSYQMTRERM',
               'O00472-205-284': 'PYRDRVIHLLALKAYKKPELLARLQKDGVNQKDKNSLGAILQQVANLNSKDLSYTLKDYVFKELQRDWPGYSEIDRRSLE',
               'O00472-522-640': 'AIELPDYLIKYIAIVSYEQRQNYKDDFNAEYDEYRALHARMETVARRFIKLDAQRKRLSPGSKEYQNVHEEVLQEYQKIKQSSPNYHEEKYRCEYLHNKLAHIKRLIGEFDQQQAESWS',
               'O00499-1-272': 'MAEMGSKGVTAGKIASNVQKKLTRAQEKVLQKLGKADETKDEQFEQCVQNFNKQLTEGTRLQKDLRTYLASVKAMHEASKKLNECLQEVYEPDWPGRDEANKIAENNDLLWMDYHQKLVDQALLTMDTYLGQFPDIKSRIAKRGRKLVDYDSARHHYESLQTAKKKDEAKIAKPVSLLEKAAPQWCQGKLQAHLVAQTNLLRNQAEEELIKAQKVFEEMNVDLQEELPSLWNSRVGFYVNTFQSIAGLEENFHKEMSKLNQNLNDVLVGLEK',
               'O00499-526-593': 'QAQHDYTATDTDELQLKAGDVVLVIPFQNPEEQDEGWLMGVKESDWNQHKELEKCRGVFPENFTERVP',
               'O00629-62-488': 'IDGDYRVQNTSLEAIVQNASSDNQGIQLSAVQAARKLLSSDRNPPIDDLIKSGILPILVHCLERDDNPSLQFEAAWALTNIASGTSEQTQAVVQSNAVPLFLRLLHSPHQNVCEQAVWALGNIIGDGPQCRDYVISLGVVKPLLSFISPSIPITFLRNVTWVMVNLCRHKDPPPPMETIQEILPALCVLIHHTDVNILVDTVWALSYLTDAGNEQIQMVIDSGIVPHLVPLLSHQEVKVQTAALRAVGNIVTGTDEQTQVVLNCDALSHFPALLTHPKEKINKEAVWFLSNITAGNQQQVQAVIDANLVPMIIHLLDKGDFGTQKEAAWAISNLTISGRKDQVAYLIQQNVIPPFCNLLTVKDAQVVQVVLDGLSNILKMAEDEAETIGNLIEECGGLEKIEQLQNHENEDIYKLAYEIIDQFFSSD',
               'O00712-1-185': 'MMYSPICLTQDEFHPFIEALLPHVRAIAYTWFNLQARKRKYFKKHEKRMSKDEERAVKDELLSEKPEIKQKWASRLLAKLRKDIRQEYREDFVLTVTGKKHPCCVLSNPDQKGKIRRIDCLRQADKVWRLDLVMVILFKGIPLESTDGERLMKSPHCTNPALCVQPHHITVSVKELDLFLAYYVQ',
               'O00716-179-350': 'YDTSLGLLTKKFIQLLSQSPDGVLDLNKAAEVLKVQKRRIYDITNVLEGIHLIKKKSKNNVQWMGCSLSEDGGMLAQCQGLSKEVTELSQEEKKLDELIQSCTLDLKLLTEDSENQRLAYVTYQDIRKISGLKDQTVIVVKAPPETRLEVPDSIESLQIHLASTQGPIEVYL',
               'O14786-25-808': 'DKCGDTIKIESPGYLTSPGYPHSYHPSEKCEWLIQAPDPYQRIMINFNPHFDLEDRDCKYDYVEVFDGENENGHFRGKFCGKIAPPPVVSSGPFLFIKFVSDYETHGAGFSIRYEIFKRGPECSQNYTTPSGVIKSPGFPEKYPNSLECTYIVFVPKMSEIILEFESFDLEPDSNPPGGMFCRYDRLEIWDGFPDVGPHIGRYCGQKTPGRIRSSSGILSMVFYTDSAIAKEGFSANYSVLQSSVSEDFKCMEALGMESGEIHSDQITASSQYSTNWSAERSRLNYPENGWTPGEDSYREWIQVDLGLLRFVTAVGTQGAISKETKKKYYVKTYKIDVSSNGEDWITIKEGNKPVLFQGNTNPTDVVVAVFPKPLITRFVRIKPATWETGISMRFEVYGCKITDYPCSGMLGMVSGLISDSQITSSNQGDRNWMPENIRLVTSRSGWALPPAPHSYINEWLQIDLGEEKIVRGIIIQGGKHRENKVFMRKFKIGYSNNGSDWKMIMDDSKRKAKSFEGNNNYDTPELRTFPALSTRFIRIYPERATHGGLGLRMELLGCEVEAPTAGPTTPNGNLVDECDDDQANCHSGTGDDFQLTGGTTVLATEKPTVIDSTIQSEFPTYGFNCEFGWGSHKTFCHWEHDNHVQLKWSVLTSKTGPIQDHTGDGNFIYSQADENQKGKVARLVSPVVYSQNSAHCMTFWYHMSGSHVGTLRVKLRYQKPEEYDQLVWMAIGHQGDHWKEGRVLLHKSLKLYQVIFEGEIGKGNLGGIAVDDISINNHISQED',
               'O14786-856-886': 'PILITIIAMSALGVLLGAVCGVVLYCACWHN',
               'Q9UJX3-1-110': 'MDPGDAAILESSLRILYRLFESVLPPLPAALQSRMNVIDHVRDMAAAGLHSNVRLLSSLLLTMSNNNPELFSPPQKYQLLVYHADSLFHDKEYRNAVSKYTMALQQKKAL',
               'Q9UJX3-133-539': 'PSEIEVKYKMAECYTMLKQDKDAIAILDGIPSRQRTPKINMMLANLYKKAGQERPSVTSYKEVLRQCPLALDAILGLLSLSVKGAEVASMTMNVIQTVPNLDWLSVWIKAYAFVHTGDNSRAISTICSLEKKSLLRDNVDLLGSLADLYFRAGDNKNSVLKFEQAQMLDPYLIKGMDVYGYLLAREGRLEDVENLGCRLFNISDQHAEPWVVSGCHSFYSKRYSRALYLGAKAIQLNSNSVQALLLKGAALRNMGRVQEAIIHFREAIRLAPCRLDCYEGLIECYLASNSIREAMVMANNVYKTLGANAQTLTLLATVCLEDPVTQEKAKTLLDKALTQRPDYIKAVVKKAELLSREQKYEDGIALLRNALANQSDCVLHRILGDFLVAVNEYQEAMDQYSIALSLD'}

mean_disorder = {'O00401': 0.5794089119018362,
                 'O00470': 0.6271861539587777,
                 'O00472': 0.5044201565435514,
                 'O00499': 0.4915315341407592,
                 'O00629': 0.14765239893593968,
                 'O00712': 0.5761809532680283,
                 'O00716': 0.5922240854771708,
                 'O14786': 0.19840866698679652,
                 'Q9UJX3': 0.15904156943216002}



def test_disorder_annotation():

    # precomputed mean disorder we can compare against


    # build a proteome
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disorder_track(P)

    for protein in P:
        assert mean_disorder[protein.unique_ID] == np.mean(protein.track('disorder').values)

    ## Test if we can change the name
    # build a proteome        
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))        
    meta_api.annotate_proteome_with_disorder_track(P, name='test')

    for protein in P:
        assert mean_disorder[protein.unique_ID] == np.mean(protein.track('test').values)

    ## Test if we raise an error if we try and overwrite
    # build a proteome        
    with pytest.raises(ProteinException):
        meta_api.annotate_proteome_with_disorder_track(P, name='test')


    ## test we can override this with safe=False
    meta_api.annotate_proteome_with_disorder_track(P, name='test', safe=False)


    ## test adding a random int for gpuid does not break things (should only be used
    # if that device is available)
    meta_api.annotate_proteome_with_disorder_track(P, name='test', safe=False, gpuid=3)
    



def test_idr_annotation():
    # Tests for the  annotate_proteome_with_disordered_domains() function in
    # the metapredict API
    #
    #

    # build a proteome
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disordered_domains(P)
    
    for d in P.domains:
        assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

    ## check we can change the name of IDRs and recall them
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disordered_domains(P, name='TEST')
    
    for d in P.domains:
        if d.domain_type == 'TEST':
            assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence
        
    ## check we can annotate folded domains as well 
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disordered_domains(P, annotate_folded_domains=True)
    
    for d in P.domains:
        if d.domain_type == 'IDR':
            assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence



    for d in P.domains:
        if d.domain_type == 'IDR':
            pass
        else:
            assert folded_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

    ## check we can rename annotate folded domains 
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disordered_domains(P, annotate_folded_domains=True, folded_domain_name='YES')
    
    for d in P.domains:
        if d.domain_type == 'IDR':
            assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

    for d in P.domains:
        if d.domain_type == 'YES':
            assert folded_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

def test_disorder_and_idr_annotation():
    # Tests for the  annotate_proteome_with_disorder_tracks_and_disordered_domains() function in
    # the metapredict API
    #
    #

    ## Test first basic annotation
    ## ................................................
    # 
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disorder_tracks_and_disordered_domains(P)


    ## check IDRs
    for d in P.domains:
        assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

    ## test disorder values
    for protein in P:
        assert mean_disorder[protein.unique_ID] == np.mean(protein.track('disorder').values)

        

    ## Test we can change names of IDRs and disorder tracks
    ## ................................................
    # 
    ## check we can change the name of IDRs and recall them
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disorder_tracks_and_disordered_domains(P, track_name='track-TEST', domain_name='domain-TEST')    
    for d in P.domains:
        if d.domain_type == 'domain-TEST':
            assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence

        assert mean_disorder[d.protein.unique_ID] == np.mean(d.protein.track('track-TEST').values)
        
    ## check we can annotate folded domains as well and we can change their name
    test_data_dir = shephard.get_data('test_data')
    P = uniprot.uniprot_fasta_to_proteome('%s/%s' % (test_data_dir, 'testset_1.fasta'))

    meta_api.annotate_proteome_with_disorder_tracks_and_disordered_domains(P, annotate_folded_domains=True, folded_domain_name='FD-test')
    
    for d in P.domains:
        if d.domain_type == 'IDR':
            assert IDR_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence
        elif d.domain_type == 'FD-test':
            assert folded_seqs[f"{d.protein.unique_ID}-{d.start}-{d.end}"] == d.sequence
