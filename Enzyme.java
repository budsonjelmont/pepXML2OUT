package pepXML2OUT_jmb;

import java.util.Vector;

public class Enzyme {
    String name;
    Vector<String> cleavage=new Vector();
    
    String getCleavage(){
        String o="";
        for(int i=0;i<cleavage.size();i++){
            o+=cleavage.get(i);
        }
        return o;
    }
}
