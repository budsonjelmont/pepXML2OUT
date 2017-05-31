package pepXML2OUT_jmb;

public class SearchHeader {
    
    // file-specific
    String filename;
    String exp_MH, charge;
    
    // search-wide
    String searchid;
    
    Modifications mods; //check terminal mod
    String db, num_of_AAs, num_of_proteins;
    String enzyme, missed_cleavage;
    String frag_tol, par_tol;
    String meta, userid, license, searchEngine;
    String display_hits, rules;
    String mass_type;
}
