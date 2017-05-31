package pepXML2OUT_jmb;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Vector;

public class Modifications {

    Vector<String[]> vModsCollection=new Vector();  // AA, delta, tag, mass 1...n
    Vector<String[]> fModsCollection=new Vector();  // AA, mass
    private String tags="*#%^@$[]~&+=|/";
    private int tagidx=1;  //* reserved for phosphorylation
    NumberFormat formatter = new DecimalFormat("#.00");
    
    String getSymbolByMass(String mass){
        System.out.println("looking for "+mass);
        for(int i=0;i<vModsCollection.size();i++){
            String[] s=vModsCollection.get(i);
            for(int j=3;j<s.length;j++){
                System.out.println("   ..."+s[j]);                
                if(formatter.format(Double.parseDouble(mass)).equals(formatter.format(Double.parseDouble(s[j])))){
                    System.out.println("matched "+s[2]); 
                    return s[2];
                }
            }
        }
        return "";
    }

    String assignTag(double dmass){
        if(Math.abs(dmass-80.00)<0.1) { //assume phosphorylation
            return String.valueOf(tags.charAt(0));
        }else
            return String.valueOf(tags.charAt(tagidx++));
    }
    
    void addVarMod(String aa, String dmass, String mass){
        for(int i=0;i<vModsCollection.size();i++){
            if(dmass.equals(vModsCollection.get(i)[1]) && aa.length()==1 && !vModsCollection.get(i)[0].equals("ct") && !vModsCollection.get(i)[0].equals("nt")){
                String[] s=vModsCollection.get(i);
                s[0]+=aa;
                String[] s1=new String[s.length+1];
                System.arraycopy(s, 0, s1, 0, s.length);
                s1[s1.length-1]=mass;
                vModsCollection.setElementAt(s1, i);
                return;
            }
        }
        String[] s=new String[4];
        s[0]=aa;
        s[1]=dmass;
        s[2]=assignTag(Double.parseDouble(dmass));
        s[3]=mass;
        vModsCollection.add(s);
    }
    
    void addFixedMod(String aa, String mass){
        String[] s=new String[2];
        s[0]=aa;
        s[1]=mass;
        fModsCollection.add(s);
    }
   

    String generateModString(){
        String modstring="";
        
        for(int i=0;i<vModsCollection.size();i++){
            modstring+="("+vModsCollection.get(i)[0]+vModsCollection.get(i)[2]+" "+vModsCollection.get(i)[1]+") ";
        }
        for(int i=0;i<fModsCollection.size();i++){
            modstring+=fModsCollection.get(i)[0]+"="+fModsCollection.get(i)[1]+" ";
        }
        return modstring;
    }
    
}
