package pepXML2OUT_jmb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import org.xml.sax.ContentHandler;
import java.util.ArrayList;
import java.util.Vector;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class DocHandlers extends DefaultHandler {

    private StringBuffer buffer = new StringBuffer();
    
    Vector<Enzyme> enzymes=new Vector();
    Vector<String> cutsites=new Vector();
    Vector<SearchHeader> searchHeader=new Vector();
    Vector<String[]> fixedMods = new Vector();
    Vector<String[]> variableMods = new Vector();
    Vector<String[]> mods = new Vector();
    
    SearchHeader sh = new SearchHeader();
    
    ArrayList<PeptideHit> peptideHit = new ArrayList<PeptideHit>();
    
    Enzyme e;
    
    Modifications m = new Modifications();
    
    private String cut, enzymeName, searchheaderID, searcheng, precursor_mass_type, fragment_mass_type,
            local_path, residue_size, size_in_db, max_intern_cleavage, enzyme, massdiff, mass, aa,
            term, prot_term, name, value, s_license, s_meta,s_userid, s_rules, s_display_hits, s_frag_tol,
            s_par_tol, spectrum, pn_mass, charge, current_search_id, index, hit_rank, peptide, ions, 
            tot_num_ions, prev_aa, next_aa, num_matched_extra_prots, protein_id, protein_name, theo_MH,
            primary_score, identscore, homoscore, expect;
        
    private String basef="", fname="";
    
    //methods to get filenames after they are generated in the main class
    public void setBaseFile(String b){ 
        basef = b;
    }
   
    public void setFileName(String f){ 
        fname = f;
    }
    
    @Override            
    
    public void startElement(String namespaceURI, String localName,
    String qName, Attributes atts) throws SAXException {
        buffer.setLength(0);
        if (localName.equals("sample_enzyme")) {
            System.out.println("sample_enzyme");
            enzymeName=atts.getValue("name");
            cutsites.clear();
        }  
        else if (localName.equals("specificity")) {
            System.out.println("specificity");
            cut = atts.getValue("cut");
            cutsites.add(cut);
        }
        else if (localName.equals("search_summary")) {
            searchheaderID = atts.getValue("search_id");
            searcheng = atts.getValue("search_engine");
            precursor_mass_type = atts.getValue("precursor_mass_type");
            fragment_mass_type = atts.getValue("fragment_mass_type");
        }
        else if (localName.equals("search_database")) {
            local_path=atts.getValue("local_path");                        
            residue_size=atts.getValue("size_of_residues");
            size_in_db=atts.getValue("size_in_db_entries");
        }
        else if (localName.equals("enzymatic_search_constraint")){
            max_intern_cleavage=atts.getValue("max_num_internal_cleavages");
            enzyme = atts.getValue("enzyme");            
        }
        else if (localName.equals("aminoacid_modification")){
            massdiff = atts.getValue("massdiff");
            mass = atts.getValue("mass");
            aa = atts.getValue("aminoacid");
            if(atts.getValue("variable").equals("Y")){
                variableMods.add(new String[]{aa, massdiff, mass});
            }    
            else if(atts.getValue("variable").equals("N")){
                fixedMods.add(new String[]{aa, mass});
                }
            }
        else if (localName.equals("terminal_modification")){
            massdiff = atts.getValue("massdiff");
            mass = atts.getValue("mass");
            term = atts.getValue("terminus");
            prot_term = atts.getValue("protein_terminus");
            if(atts.getValue("variable").equals("Y")){
                variableMods.add(new String[]{term.toLowerCase()+"t", massdiff, mass});
            }    
            else if(atts.getValue("variable").equals("N")){
                fixedMods.add(new String[]{prot_term.isEmpty()?
                        "+"+term.toUpperCase()+"term-pep":
                        "+"+term.toUpperCase()+"term-prot", 
                        mass});
                }
            }
        else if (localName.equals("parameter")){
            name = atts.getValue("name");
            value = atts.getValue("value");
            if(name.equalsIgnoreCase("LICENSE"))
                    s_license=value;
                if(name.equalsIgnoreCase("COM"))
                    s_meta=value;
                if(name.equalsIgnoreCase("USERID"))
                    s_userid=value;
                if(name.equalsIgnoreCase("RULES"))
                    s_rules=value;
                if(name.equalsIgnoreCase("REPORT"))
                    s_display_hits=value;
                if(name.equalsIgnoreCase("ITOL"))
                    s_frag_tol=value;
                if(name.equalsIgnoreCase("TOL"))
                    s_par_tol=value;
        }
        else if(localName.equals("spectrum_query")){ //if a spectrum query has no search results, then it will have no child nodes
            //empty the peptide hits arraylist to prepare for new results 
            peptideHit.clear();
            spectrum = atts.getValue("spectrum");
            System.out.println("Results for spectrum: " + spectrum);
            pn_mass = atts.getValue("precursor_neutral_mass");
            charge = atts.getValue("assumed_charge");
            index = atts.getValue("index");
        }
        else if(localName.equals("search_result")){
            current_search_id = atts.getValue("search_id");
            SearchHeader s = null;
            //determine which searchcheader is appropriate for this search
              for(int k=0;k<searchHeader.size();k++){
                if(searchHeader.get(k).searchid.equals(current_search_id)){
                    s=searchHeader.get(k);
                    break;
                }
              }
              if(s==null&&searchHeader.size()==1){
                    s=searchHeader.get(0);
              }
              if(s!=null){
                String sprm=spectrum;
                s.charge="+"+charge;
              // determine filename
                if(sprm.indexOf(".dta")>=0)
                    s.filename=sprm.substring(0, sprm.lastIndexOf("."))+".out";
                else if(sprm.length()>9) //assume it is all dta name except .dta
                   s.filename=sprm+".out";
                else if(isAllNumber(sprm)) //scan number
                   s.filename=fname+"."+sprm+"."+sprm+"."+s.charge.substring(1)+".out";
                else if(sprm.indexOf("spectrumId")>=0){
                    String sn=getNextNum(sprm,sprm.indexOf("spectrumId")+"spectrumId".length());
                    s.filename=fname+"."+sn+"."+sn+"."+s.charge.substring(1)+".out";
                }else if(sprm.indexOf("cqNumber")>=0){
                    String sn=getNextNum(sprm,sprm.indexOf("cqNumber")+"cqNumber".length());
                    s.filename=fname+"."+sn+"."+sn+"."+s.charge.substring(1)+".out";
                }else
                    s.filename=fname+"."+index+"."+index+"."+s.charge.substring(1)+".out";
                s.exp_MH=String.valueOf(Double.parseDouble(pn_mass)+1.007274);
            }
            //assign the searcheader specified in this method to the searchheader sh that is shared in this class 
            sh = s;  
        }      
        else if(localName.equals("search_hit")){
	        	
            hit_rank=atts.getValue("hit_rank");
            peptide=atts.getValue("peptide");
            ions=atts.getValue("num_matched_ions")+"/"+
            (atts.getValue("tot_num_ions")==null||atts.getValue("tot_num_ions").isEmpty()
                ? 2*(peptide.length()-1) 
                : atts.getValue("tot_num_ions"));
            prev_aa= atts.getValue("peptide_prev_aa")==null||atts.getValue("peptide_prev_aa").isEmpty()?"-":atts.getValue("peptide_prev_aa");
            next_aa= atts.getValue("peptide_next_aa")==null||atts.getValue("peptide_next_aa").isEmpty()?"-":atts.getValue("peptide_next_aa");
            num_matched_extra_prots=Integer.parseInt(atts.getValue("num_tot_proteins"))>1?"+"+(Integer.parseInt(atts.getValue("num_tot_proteins"))-1):"";
            protein_id=atts.getValue("protein");
            protein_id+="|";
            String prefix="";
            if(protein_id.startsWith("#")){
                int ps=protein_id.lastIndexOf("#");
                prefix=protein_id.substring(0,ps+1);
                protein_id=protein_id.substring(ps);
            }
            if(protein_id.startsWith("IPI") && !protein_id.substring(3,4).equals(":")){
                protein_id="IPI:"+protein_id.substring(3);
            }else if(protein_id.startsWith("gi") && !protein_id.substring(2,3).equals("|")){
                protein_id="gi|"+protein_id.substring(2);
            }
            protein_id=prefix+protein_id;
            protein_name=atts.getValue("protein_descr")==null||atts.getValue("protein_descr").isEmpty()?
                protein_id.substring(0,protein_id.length()-1):
                atts.getValue("protein_descr");
            NumberFormat formatter = new DecimalFormat("#.00000");
            theo_MH=String.valueOf(formatter.format(Double.parseDouble(atts.getValue("calc_neutral_pep_mass"))+1.007274));
        }
        else if(localName.equals("search_score")){
            if(searcheng.equalsIgnoreCase("MASCOT")){
                if(atts.getValue("name").equalsIgnoreCase("ionscore"))
                    primary_score=atts.getValue("value");
                if(atts.getValue("name").equalsIgnoreCase("identityscore"))
                    identscore=atts.getValue("value");
                if(atts.getValue("name").equalsIgnoreCase("homologyscore"))
                    homoscore=atts.getValue("value");
                if(atts.getValue("name").equalsIgnoreCase("expect")){
                    expect=atts.getValue("value");
                    double d=Double.parseDouble(expect);
                    d=-1*Math.log10(d);
                    d=Math.round(d*1000)/1000d;
                    expect=String.valueOf(d);
            }
            } else if(searcheng.equalsIgnoreCase("OMSSA")||true){
                if(atts.getValue("name").equalsIgnoreCase("pvalue")){
                    expect=atts.getValue("value");
                    double d=Double.parseDouble(expect);
                    d=-1*Math.log10(d);
                    d=Math.round(d*1000)/1000d;
                    expect=String.valueOf(d);
                }
                if(atts.getValue("name").equalsIgnoreCase("expect")){
                    primary_score=atts.getValue("value");
                    double d=Double.parseDouble(primary_score);
                    d=-1*Math.log10(d);
                    d=Math.round(d*1000)/1000d;
                    primary_score=String.valueOf(d);
                    //if(expect.indexOf("e")>0){
                    //    String[] tmp=expect.split("e");
                    //    tmp[1]=tmp[1].substring(1,tmp[1].length());
                    //    expect=String.valueOf(Double.parseDouble(tmp[0])*(10^(Integer.parseInt(tmp[1]))));
                   // }
                }
            }
        }
      else if (localName.equals("modification_info")){
            //mods.clear(); 130911
            if(atts.getValue("mod_cterm_mass")!=null && !atts.getValue("mod_cterm_mass").isEmpty()){
                String[] mod=new String[2];
                mod[0]=String.valueOf(peptide.length());
                mod[1]=sh.mods.getSymbolByMass(atts.getValue("mod_cterm_mass"));
                mods.add(mod);
            }
            if(atts.getValue("mod_nterm_mass")!=null && !atts.getValue("mod_nterm_mass").isEmpty()){
                String[] mod=new String[2];
                mod[0]=String.valueOf(1);
                mod[1]=sh.mods.getSymbolByMass(atts.getValue("mod_nterm_mass"));
                mods.add(mod);
            }
      }
       else if (localName.equals("mod_aminoacid_mass")){
                String[] mod=new String[2];
                mod[0]=atts.getValue("position");
                mod[1]=sh.mods.getSymbolByMass(atts.getValue("mass"));
                mods.add(mod);
       }        
    }
 @Override
        
 public void endElement(String uri, String localName, String qName)throws SAXException {      
  if (localName.equals("sample_enzyme")) {
        //System.out.println("temp val:");  
        //System.out.println(buffer.toString());
        Enzyme e=new Enzyme();
        e.name = enzymeName;
        for(int i = 0; i < cutsites.size(); i++){
            System.out.println("cut sites:");
            System.out.println(cutsites.elementAt(i));
            e.cleavage.add(cutsites.elementAt(i));
        }
        enzymes.add(e);
  }

  else if (localName.equals("search_summary")) {
      SearchHeader s=new SearchHeader();
      s.searchid=searchheaderID;                    
      s.searchEngine= searcheng==null||searcheng.isEmpty()?"Unknown Search Engine":searcheng;
      s.mass_type= precursor_mass_type==null ||precursor_mass_type.isEmpty() ? "UNKNOWN": (precursor_mass_type.startsWith("mono")?"MONO":"AVE");
      s.mass_type=s.mass_type+"/"+ (fragment_mass_type==null ||fragment_mass_type.isEmpty() ? "UNKNOWN": (fragment_mass_type.startsWith("mono")?"MONO":"AVE"));
      s.db= local_path==null||local_path.isEmpty()?"UNKNOWN":local_path.replaceAll("/", "\\\\");
      s.num_of_AAs = residue_size==null||residue_size.isEmpty()?"######":residue_size;
      s.num_of_proteins= size_in_db==null||size_in_db.isEmpty()?"###":size_in_db;
      s.missed_cleavage= max_intern_cleavage==null||max_intern_cleavage.isEmpty()?"N/A":max_intern_cleavage;
      s.enzyme= enzyme==null||enzyme.isEmpty()?"UNKNOWN":enzyme;
      for(int k=0;k<enzymes.size();k++){
         if(enzyme.equals(enzymes.get(k).name)){
            s.enzyme=enzymes.get(k).name+"("+enzymes.get(k).getCleavage()+")";
            break;
         }
      }
      s.license = s_license;
      s.meta = s_meta;
      s.userid = s_userid;
      s.rules = s_rules;
      s.display_hits = s_display_hits;
      s.frag_tol = s_frag_tol;
      s.par_tol = s_par_tol;
      s.license= s.license==null || s.license.isEmpty() ? "LICENSE: NO INFO" : s.license;
      s.meta= s.meta==null || s.meta.isEmpty() ? "SEARCH SUBMITTED" : s.meta;
      s.userid= s.userid==null || s.userid.isEmpty() ? "UNKNOWN" : s.userid;
      s.rules= s.rules==null || s.rules.isEmpty() ? "N/A" : s.rules;
      s.display_hits= s.display_hits==null || s.display_hits.isEmpty() ? "DEFAULT" : s.display_hits;
      s.frag_tol= s.frag_tol==null || s.frag_tol.isEmpty() ? "ND" : (s.frag_tol.substring(0,1).equals(".")?"0"+s.frag_tol:s.frag_tol);
      s.par_tol= s.par_tol==null || s.par_tol.isEmpty() ? "ND" :(s.par_tol.substring(0,1).equals(".")?"0"+s.par_tol:s.par_tol);
      s.mods=new Modifications();
      //loop through temp vectors to get all the variable and fixed modifications we captured
      System.out.println("Adding variable mods to search header");
      for(int i = 0; i < variableMods.size(); i++){
          s.mods.addVarMod(variableMods.get(i)[0], variableMods.get(i)[1], variableMods.get(i)[2]);
          System.out.println("amino acid: " + variableMods.get(i)[0]);
          System.out.println("mass difference: " + variableMods.get(i)[1]);
          System.out.println("mass: " + variableMods.get(i)[2]);
      }
      System.out.println("Adding fixed mods to search header");
      for(int i = 0; i < fixedMods.size(); i++){
          s.mods.addFixedMod(fixedMods.get(i)[0], fixedMods.get(i)[1]);
          System.out.println("amino acid: " + fixedMods.get(i)[0]);
          System.out.println("mass : " + fixedMods.get(i)[1]);
      }
      //System.out.println(s.license);
      searchHeader.add(s);
      //clear out temporary storage vectors so they are ready to accept new parameters if there's another search header
      variableMods.clear();
      fixedMods.clear();
  }

  else if (localName.equals("search_hit")) {

    PeptideHit hit=new PeptideHit();
    hit.rank = hit_rank;
    hit.peptide = peptide;
    hit.ions = ions;
    hit.prevAA = prev_aa;
    hit.nextAA = next_aa;
    hit.num_matched_extra_prots = num_matched_extra_prots;
    hit.protein_id = protein_id;
    hit.protein_name = protein_name;
    hit.theo_MH=theo_MH;
    
    if(searcheng.equalsIgnoreCase("MASCOT")){
        hit.primary_score = primary_score;
        hit.identscore = identscore;
        hit.homoscore = homoscore;;
        hit.expect = expect;
    } else if(searcheng.equalsIgnoreCase("OMSSA")||true){
        hit.primary_score = primary_score;
        hit.expect = expect;
    }
    if(mods.size()>0){   
            hit.tagged_sequence=hit.prevAA+"."+hit.getSeqWithMods(mods)+"."+hit.nextAA;
        }else{
            hit.tagged_sequence=hit.prevAA+"."+hit.peptide+"."+hit.nextAA;
        }
    peptideHit.add(hit);
    	
    mods.clear();
  }
  else if (localName.equals("spectrum_query")) {
//      SearchHeader s = null;
//      for(int k=0;k<searchHeader.size();k++){
//        if(searchHeader.get(k).searchid.equals(current_search_id)){
//            s=searchHeader.get(k);
//            break;
//        }
//      }
//      if(s==null&&searchHeader.size()==1){
//            s=searchHeader.get(0);
//      }
//      if(s!=null){
//        String sprm=spectrum; //Already read in spectrum
//        s.charge="+"+charge;
//      // determine filename
//        if(sprm.indexOf(".dta")>=0)
//            s.filename=sprm.substring(0, sprm.lastIndexOf("."))+".out";
//        else if(sprm.length()>9) //assume it is all dta name except .dta
//           s.filename=sprm+".out";
//        else if(isAllNumber(sprm)) //scan number
//           s.filename=fname+"."+sprm+"."+sprm+"."+s.charge.substring(1)+".out";
//        else if(sprm.indexOf("spectrumId")>=0){
//            String sn=getNextNum(sprm,sprm.indexOf("spectrumId")+"spectrumId".length());
//            s.filename=fname+"."+sn+"."+sn+"."+s.charge.substring(1)+".out";
//        }else if(sprm.indexOf("cqNumber")>=0){
//            String sn=getNextNum(sprm,sprm.indexOf("cqNumber")+"cqNumber".length());
//            s.filename=fname+"."+sn+"."+sn+"."+s.charge.substring(1)+".out";
//        }else
//            s.filename=fname+"."+index+"."+index+"."+s.charge.substring(1)+".out";
//        s.exp_MH=String.valueOf(Double.parseDouble(pn_mass))+1.007274;
//      
//        //now that we know which searchheader we will be using, call the getSymbolByMass method
//        //on the arguments read in from the mod_aminoacid tags
//        for(int i = 0 ; i < mods.size(); i++){
//            mods.get(i)[1] = s.mods.getSymbolByMass(mods.get(i)[1]);
//        }
//        if(mods.size())>0){   
//            peptideHit[k].tagged_sequence=peptideHit[k].prevAA+"."+peptideHit[k].getSeqWithMods(mods)+"."+peptideHit[k].nextAA;
//        }else{
//            peptideHit[k].tagged_sequence=peptideHit[k].prevAA+"."+peptideHit[k].peptide+"."+peptideHit[k].nextAA;
//        }
//        writeToFile(sh, peptideHit);
//      }
  }
   else if(localName.equals("search_result")){
       writeToFile(sh, peptideHit);
   }
}

    private String tmpValue;
    
    @Override
    public void characters(char[] ch, int length, int start) throws SAXException {
        buffer.append(ch, start, length);
        //System.out.println(start);
        //System.out.println(length);
        //System.out.println(">"+buffer.toString());
    }   
       void writeToFile(SearchHeader s, ArrayList<PeptideHit> ph){
        //System.out.println(basef+s.filename);
        try{
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(basef+s.filename)));
            out.write(System.getProperty("line.separator"));
            out.write(" "+s.filename+System.getProperty("line.separator"));
            out.write(" "+s.searchEngine+System.getProperty("line.separator"));
            if(s.searchEngine.equalsIgnoreCase("MASCOT"))
                out.write(" Field redefined: XCorr->Ion Score(Mowse Score); Sp->-log10(expect value)"+System.getProperty("line.separator"));
            else
                out.write(" Field redefined: XCorr->-log10(expect value); Sp->-log10(p-value)"+System.getProperty("line.separator"));
            out.write(" "+s.license+System.getProperty("line.separator"));
            out.write(" "+s.meta+"; User: "+s.userid+System.getProperty("line.separator"));
            out.write(" (M+H)+ mass = "+s.exp_MH+" ~ "+s.par_tol+" ("+s.charge+"), fragment tol = "+s.frag_tol+" , "+s.mass_type+System.getProperty("line.separator"));
            out.write(" NADA"+System.getProperty("line.separator"));
            out.write(" # amino acids = "+s.num_of_AAs+", # proteins = "+s.num_of_proteins+", "+s.db+System.getProperty("line.separator"));
            out.write(" Rules: "+s.rules+System.getProperty("line.separator"));
            out.write(" display top "+s.display_hits+"/5, ion % = 0.0, CODE = N/A"+System.getProperty("line.separator"));
            out.write(" "+s.mods.generateModString()+" Enzyme:"+s.enzyme+" ("+s.missed_cleavage+")"+System.getProperty("line.separator"));
            out.write(System.getProperty("line.separator"));

            int maxRef=0, maxExtraPro=0;
            double maxIonScore=0.0;
            for(int i=0;i<ph.size() && i<10;i++){
                if(ph.get(i).protein_id.length()>maxRef)
                    maxRef=ph.get(i).protein_id.length();
                if(Double.parseDouble(ph.get(i).primary_score)>maxIonScore)
                    maxIonScore=Double.parseDouble(ph.get(i).primary_score);
                if(ph.get(i).num_matched_extra_prots.length()>maxExtraPro)
                    maxExtraPro=ph.get(i).num_matched_extra_prots.length();                
            }
            int fill="Reference".length()>maxRef?0:maxRef-"Reference".length();
            fill+=maxExtraPro>0?maxExtraPro+3:2;
            out.write("  #   Rank/Sp      Id#     (M+H)+    deltCn   XCorr    Sp    Ions   Reference");
            out.write(Space(fill)+"Peptide"+System.getProperty("line.separator"));
            out.write(" ---  --------  --------  --------   ------  ------   -----  -----  ---------");
            out.write(Space(fill)+"-------"+System.getProperty("line.separator"));
            
            NumberFormat formatter = new DecimalFormat("0.0000");
            for(int i=0;i<ph.size() && i<10;i++){
                out.write(" ");
                out.write(strFormat((i+1)+".",3,false));
                out.write("  ");
                out.write(strFormat(ph.get(i).rank+" / ## ",8,false));
                out.write("  ");
                out.write(strFormat("N/A",8,false));
                out.write(" ");
                out.write(strFormat(ph.get(i).theo_MH,10,false));
                out.write("  ");
                out.write(strFormat(formatter.format((maxIonScore-Double.parseDouble(ph.get(i).primary_score))/maxIonScore),6,false));
                out.write("  ");
                out.write(strFormat(ph.get(i).primary_score,6,false));
                out.write("  ");
                out.write(strFormat(ph.get(i).expect,6,false));
                out.write("  ");
                out.write(strFormat(ph.get(i).ions,5,false));
                out.write("  ");
                out.write(ph.get(i).protein_id);
                out.write(Space(1+("Reference".length()>maxRef?"Reference".length()-ph.get(i).protein_id.length():maxRef-ph.get(i).protein_id.length())));
                out.write(strFormat(ph.get(i).num_matched_extra_prots,maxExtraPro+2,true));
                out.write(ph.get(i).tagged_sequence);
                out.write(System.getProperty("line.separator"));
            }
            out.write(System.getProperty("line.separator"));
            for(int i=0;i<ph.size() && i<5;i++){
                out.write(" ");
                out.write(strFormat((i+1)+".",3,false));
                out.write("  ");
                out.write(ph.get(i).protein_id);
                out.write(" ");
                out.write(ph.get(i).protein_name);
                out.write(System.getProperty("line.separator"));
            }
            
            out.close();
        }catch(Exception e){
            e.printStackTrace();
        }
    }
    
    String strFormat(String s, int width, boolean leftAlign){
        if(s==null)
            s="N/A";
        if(s.length()<width){
            if(leftAlign)
                s=s+Space(width-s.length());
            else
                s=Space(width-s.length())+s;
        }else if(s.length()>width){
            s=s.substring(0,width);
        }
        return s;
    }
    
    String Space(int num){
        String o="";
        for(int i=0;i<num;i++){
            o+=" ";
        }
        return o;
    }
    
    boolean isAllNumber(String input){
        boolean valid=false;
        try{
            int x=Integer.parseInt(input);
            valid=true;
        }catch(Exception e){
            valid=false;
        }
        return valid;
    }
    
    String getNextNum(String str, int idx){
        StringBuilder sb=new StringBuilder();
        boolean beginAdd=false;
        for(int i=idx;i<str.length();i++){
            if(isAllNumber(str.substring(i, i+1))){
                if(!beginAdd)
                    beginAdd=true;
                sb.append(str.substring(i, i+1));
            }else{
                if(beginAdd)
                    break;
            }
        }
        return sb.toString();
    }
  
    void reportError(String msg){
        System.out.println("Error: " + msg);
        try{
            //BufferedWriter out = new BufferedWriter(new FileWriter(new File(AppFolder+"mascotdownloaderoutput.txt")));
            //out.write("Error: " + msg);
            //out.close();
        }catch(Exception e){
            
        }
    }
   
private static String convertToFileURL(String filename) {
        String path = new File(filename).getAbsolutePath();
        if (File.separatorChar != '/') {
            path = path.replace(File.separatorChar, '/');
        }

        if (!path.startsWith("/")) {
            path = "/" + path;
        }
        return "file:" + path;
    }   
}