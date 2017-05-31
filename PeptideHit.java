package pepXML2OUT_jmb;

import java.util.Vector;

public class PeptideHit {
    String rank;
    String peptide, tagged_sequence;
    String protein_id, protein_name, num_matched_extra_prots;
    String prevAA, nextAA;
    String theo_MH;
    String ions;
    String primary_score, identscore, homoscore, expect; // mascot type
    /**
    String xcorr, sp, deltaCn;                      // sequest type
    */
    String getSeqWithMods(Vector<String[]> modstr){
        if(modstr.size()==0)
            return peptide;
        else{
            // modstr: 0-pos, 1-symbol
            int len=peptide.length();
            String[] s=new String[len];
            for(int i=0;i<len;i++){
                s[i]=String.valueOf(peptide.charAt(i));
            }
            for(int i=0;i<modstr.size();i++){
                int pos=Integer.parseInt(modstr.get(i)[0]);
                if(pos<=len && pos>0){
                    s[pos-1]+=modstr.get(i)[1];
                }
            }
            String output="";
            for(int i=0;i<len;i++){
                output+=s[i];
            }
            return output;
        }
    }
}
