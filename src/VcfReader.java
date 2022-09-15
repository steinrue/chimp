import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.CloseableIterator ;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader ;

public class VcfReader {
	
	///////////////////////////////////////////////////////////////
	// MAIN METHOD ////////////
	////////////////////////////////
	public static void main(String[] args) {
		Helper helper = new Helper();
		
		// TODO Auto-generated method stub
		String folder = "/Users/gautam/Desktop/PD_diag/";
		VcfReader vRead = new VcfReader(folder + "GBR_Chr2.vcf", folder + "chr2.ref", folder + "homo_sapiens_ancestor_2.fa", new int[] {120000000,150000000}, new int[] {1,10});
		//vRead.update_max_tract_size(20000);
		
		vRead.print_stream_to_file(folder + "stream.csv");
		
		//int[][] stream = vRead.getStream(500);
		

		System.out.println("\n DONE!!!");
	}

	
	
	// all of inputs needed to specify how we extract the output stream  //
	File target_vcf_file;
	String reference_file;
	String ancestral_file;
	
	int chrom_seg_lb; // left bound position of the chromosomal sequence we want to consider
	int chrom_seg_rb;  // right boundary of sequence..
	
	int num_samples; // total num of sample indices included in the interval between lb, rb
	int[] hap_indices; // the indices for the haplotypes to create the stream for start with first individual, count all their haps, continue
	
	int max_tract_size = 1000;
	
	
	
	// counters for various variant combos
	
	int sites_missing_from_ref; 		// missing as tagged in ref
	int sites_mismatch_from_anc; 		// anc doesnt match ref/alt from SNP_unit (could be missing)
	int sites_missing_ref_notVcf;	// missing in reference, but variant appears in vcf
	int non_snps; 					// variant is not SNP in vcf
	int snps_missing_from_ancestral; // valid snp, but missing in ancestral
	int snps_missing_from_vcf; 		// valid ref,snp,anc, but in vcf fully listed as missing
	int snps_missing_partial_from_vcf; // same as above but partially missing in vcf
	int snps_multiallelic; 			// more than just reference and 1 alternate allele (doesnt count missing)
	int valid_snps; 				// the SNPs that pass the filters, and have segregating sites
	int vcf_variant_count;
	int vcf_duplicate_position;
	
	VcfReader(String a, String b, String c, int[] chrom_bounds, int[] haplotypes){
		target_vcf_file = new File(a);
		reference_file = b;
		ancestral_file = c;
		
		num_samples = haplotypes.length;
		hap_indices = haplotypes;

		chrom_seg_lb = chrom_bounds[0];
		chrom_seg_rb = chrom_bounds[1]; 
		
		
		if( chrom_seg_lb > chrom_seg_rb ) {System.err.print("Invalid chromosome segment or sample index bounds"); System.exit(0);}
	}
	
	VcfReader(String a, String b, String c, int chrom_size, int[] hap_indices){		
		this(a,b,c, new int[] {1, chrom_size > 0 ? chrom_size : check_lengths(b,c)} , hap_indices);
	}


	///////////////////////
	// PRIMARY METHODS ////
	
	private int[][] getReducedStream() throws FileNotFoundException{
		double stime = System.nanoTime(); // just to time this method
		
		// setup holder arraylists
		ArrayList<Integer> tract_ss = new ArrayList<Integer>(0);
		ArrayList<Integer> tract_pos = new ArrayList<Integer>(0);
		ArrayList<Integer> tract_type = new ArrayList<Integer>(0);
		
		// initialize variant counters
		
		sites_missing_from_ref = 0; 		
		sites_missing_ref_notVcf = 0;	
		non_snps = 0; 					
		snps_missing_from_ancestral = 0; 
		snps_missing_from_vcf = 0; 		
		snps_missing_partial_from_vcf = 0; 
		snps_multiallelic = 0; 
		valid_snps = 0; 
		vcf_variant_count = 0;
		vcf_duplicate_position = 0;
		
		
		//****IMPLEMENT SOMETHING TO ENSURE THE SEG SITES ARE IN CORRECT ORDER
		
	
		// give a nice error message if the file does not exist
		if (!target_vcf_file.exists() || !target_vcf_file.canRead()) {
			throw new FileNotFoundException(target_vcf_file.getPath() + " (No such file or directory, or insufficient permissions)");
		}
		
		// get the variant-reading iterator setup
		VCFFileReader r_vcf = new VCFFileReader(target_vcf_file, false);
		CloseableIterator<VariantContext> vcf_iter = r_vcf.iterator();
		
		
		// setup fasta streams for reference and ancestral files, fastforward to chrom_seg_lb
		FastaStream ref_stream = new FastaStream(reference_file, chrom_seg_lb);
		FastaStream anc_stream = new FastaStream(ancestral_file, chrom_seg_lb);
		
		
		
		//NOW READY TO GO BASE BY BASE
		
		// collect first bases @ lb, and first variant at or after lb
		int position = chrom_seg_lb;
		char next_anc_base = anc_stream.next_char();
		char next_ref_base = ref_stream.next_char();
		VariantContext next_var = vcf_iter.next();
		while(next_var.getStart() < chrom_seg_lb) {next_var = vcf_iter.next();} // fastforward in vcf to line in our chromosomal target region
		

		int t_type=10; // 1 is segregating, 2 is missing
		int last_info_pos = position; // last position with info from VCF
		
		// PROCESS FIRST BASE, and start off first tract
		tract_pos.add(position);
		int[] base_info = this.process_position(position, next_ref_base, next_anc_base, next_var);
		if(base_info[0] != 0 ) { // classify first tract as missing if any number of missing elements
			tract_ss.add(-1);
			tract_type.add(2);
			t_type = 2;
		}
		else {	// first base is not missing as per vcf or ref file or anc file
			tract_ss.add(base_info[1]);
			tract_type.add(1);
			t_type = 1;

		}
		base_info = null;
		
		
		
		// NOW ITERATE ALONG THE REF POSITIONS AND PROCESS BASES
		for(position = chrom_seg_lb+1 ; position < chrom_seg_rb + 1 ; position++) {
			next_anc_base = anc_stream.next_char();
			next_ref_base = ref_stream.next_char();
			
			//make sure that the next variant is not behind the current position, count variants with same position
			while(next_var != null && next_var.getStart() < position) { 
				int prev_pos = next_var.getStart();
				next_var = vcf_iter.next(); vcf_variant_count++; 
				
				// CHECK IF DUPLICATE (if duplicate, treat as missing even if first line could be valid snp line.)
				if(next_var != null && next_var.getStart() == prev_pos) { // CHECK next_var is duplicate position of prev_pos
					vcf_duplicate_position++; // counter
				
					/*
					// check if latest tract is tied to duplicate entry, retroactively change that tract to missing 
					if(tract_pos.get(tract_pos.size()-1) == prev_pos){
						tract_type.set(tract_type.size()-1, 2);
						tract_ss.set(tract_ss.size()-1,-1);
					}
					
					// otherwise create a new missing tract
		
					else if(t_type == 2){} // if continuation, its already accounted for
					
					else{
						tract_pos.add(last_info_pos + 1);
						tract_type.add(2);
						tract_ss.add(-1);
					}
					*/
				}
				
				// this skips variant lines that have already appeared, we handle skipping the first of such a sequence shortly
			} 
			
			// make sure didn't reach end of file in ref or anc
			if(next_anc_base == (char) 0 || next_ref_base == (char) 0) {System.out.println("Trying to read past end of ref/anc file."); System.exit(0);}
			
			
			base_info = process_position(position, next_ref_base, next_anc_base, next_var);
			
			
			// if any elements are missing handle as missing
			if(base_info[0] != 0) { 
				if(t_type == 2) { /* just continuation of missing tract*/}
				else { // start of new missing tract, previous tract was part of ss style tract
					
					t_type = 2;
					
					// missing tag stems from reference
					if(next_ref_base == 'n' || next_ref_base == 'N') { 		
						tract_pos.add(position); 
						tract_ss.add(-1);
						tract_type.add(2);
					}
					
					// however if its missing as a result of non_process_missing, we tag the tract to the left (till we hit another site with info) as missing
					else {  // last_info_pos is most recent info position from vcf, 
						// we use most recent pos, not necessarily last_info_pos though
						// because last_info_pos gets its info from vcf, 
						// but tract_pos can have more recent info if ref-file indicates a missing segment recently
						// or if ending a missing tract started a 0-monomorphic tract recently
						
						// if previous tract is not starting from an entry of the VCF, convert it to missing as part of painting 
						if (tract_pos.get(tract_pos.size()-1) > last_info_pos) {
							tract_type.set(tract_type.size()-1, 2);
							tract_ss.set(tract_ss.size()-1, -1);
						}
						
						else { // previous tract is directly tied to vcf variant, so paint AFTER that position till current pos
							tract_pos.add( last_info_pos + 1); 
							tract_ss.add(-1);
							tract_type.add(2);
						}
						
					}
				}
				
				
			}
			
			// otherwise treat as non_missing
			else { 
				// if all non-missing are non-ancestral alleles, then convert to 0 non-ancestral
				if(base_info[1] == num_samples - base_info[0]) {base_info[1] = 0 ;}
				
				if(t_type == 2) {// switch from missing to non_missing tract type
					tract_ss.add(base_info[1]);
					tract_pos.add(position);
					tract_type.add(1);
					t_type = 1;
				}
				
				else { // previous tract is not missing, check if continuation
					if(base_info[1] == 0 ){/* continuation of no_seg_site tract */}
					else { // this is seg site, start new tract
						tract_ss.add(base_info[1]);
						tract_pos.add(position);
						tract_type.add(1);
						t_type = 1;
					}
					
				}
				
			}
			
			
			//update last_info_pos
			if(base_info[2] == 1) {last_info_pos = position;}

		}
		
		
		// CHECK TO HANDLE LAST POSITION IN SEQUENCE
		
		if(tract_pos.get(tract_ss.size()-1) != chrom_seg_rb) { // handle final tract (size 1) at last locus
			tract_pos.add(chrom_seg_rb);
			if(t_type == 1) {tract_ss.add(0); tract_type.add(1); }
			if(t_type == 2) {tract_ss.add(-1); tract_type.add(2); }
		}
		
		
		

		
	
		
		
		///////////////////////////////////////////////////
		//	PROCESS AND CREATE OUTPUT STRUCTURE
		/////////////////////////////////////////////////
		int n_tracts = tract_ss.size();
		
		// ADD IN ADDITIONAL TRACTS TO MEET MAX_TRACT_SIZE REQUIREMENT
		int next_default_position = chrom_seg_lb + max_tract_size;
		for(int i = 0 ; i < n_tracts; i++) {
			
			// check if we've reached cutting point
			if(tract_pos.get(i) > next_default_position) {
				
				tract_pos.add(i, next_default_position);
				// duplicate the previous tract type at current position
				if(tract_type.get(i-1) == 1) {
					tract_type.add(i,1);
					tract_ss.add(i,0);
				}
				else if(tract_type.get(i-1) == 2) {
					tract_type.add(i,2);
					tract_ss.add(i,-1);
				}
				n_tracts++;	
			}
	
			next_default_position = tract_pos.get(i) + max_tract_size;
		}
		
		
		
		
		// create OUTPUT ARRAYS (convert from position to tract lengths)
		n_tracts = tract_ss.size();
		int[][] ssStream = new int[3][n_tracts];
			
		for(int i = 0; i < n_tracts; i++) {
			   ssStream[0][i] = tract_ss.get(i);
			   ssStream[1][i] = tract_pos.get(i);
			   ssStream[2][i] = tract_type.get(i);
				   
			   if(ssStream[0][i] > 0) {valid_snps++;}
			  
			   if(i== n_tracts - 1) {ssStream[1][i] = 1;}
			   else { ssStream[1][i] = tract_pos.get(i+1) - tract_pos.get(i);}
				   
		}
		
		
		
		try {
			ref_stream.buff_read.close();
			anc_stream.buff_read.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		r_vcf.close();
		
		
		// Print stats and close buffers
		System.out.println("File Read: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds");
		System.out.println("Read in " + this.target_vcf_file.getName() + "  ::  " + " haplotypes={"+Arrays.toString(hap_indices)+"} over " + (chrom_seg_rb - chrom_seg_lb) + " loci :: *** " + valid_snps  + " valid SNPs ***");


		
		return ssStream;
	}
	

	int[] process_position(int pos, char ref_ch , char anc_ch, VariantContext variant) {
		
		int missing_haps = 0;			// params to return	(number of missing entries at pos)
		int non_ancestral_haps = 0;		// params to return (number of derived alleles at this pos)
		int info_at_position = 0;		// param to return (0 if no VCF entry for this pos, 1 otherwise)
		
		// if REF says MISSING, then tag as missing and move on
		if(ref_ch == 'N' || ref_ch == 'n') { 
			missing_haps = num_samples; 
			sites_missing_from_ref++;
			if(variant != null && variant.getStart() == pos) {sites_missing_ref_notVcf++;} // count as missing in reference but not in vcf
		} // ref is not missing
		
		// if NO VAR at site, assume its non segregating and move on
		else if(variant == null || variant.getStart() != pos) { 
			missing_haps = 0; non_ancestral_haps = 0;
		} 
		
		//  HAVE VALID VARIANT AT SITE
		
		else {
			info_at_position = 1;

			// make sure REF and VCF files are SYNCHED properly, else quit
			if(variant.getReference().getBaseString().charAt(0) != ref_ch ) {
				System.out.println(variant.getReference().getBaseString().charAt(0) + "!=" + ref_ch + " at position "+ pos );
				System.out.println("VCF doesnt synch with reference file"); System.exit(0);
			}
			
			// get info from variant
			SNP_unit snp = new SNP_unit(variant);
			missing_haps = snp.getMissing(); non_ancestral_haps = snp.getDerived();

			String anc_str = Character.toString(anc_ch);
			
			// make sure ANC is NONMISSING - else tag as missing
			if( isATCG(anc_ch) == false ) {
				snps_missing_from_ancestral++; missing_haps = num_samples; 
			}
			
			// PASS IF NON-SEGREGATING FOR RELEVANT SAMPLES
			else if( snp.no_alternate && !snp.part_missing && !snp.full_missing) { // effectively non-segregating
				missing_haps = 0; non_ancestral_haps = 0;
			}
			
			// IF VARIANT is NOT SNP, tag as missing and move on
			else if(variant.getType().equals(VariantContext.Type.SNP) == false){
				non_snps++; missing_haps = num_samples;
			}
			
			// if segregating, but multiallelic, tag as missing
			else if(snp.is_multiallelic) {
				snps_multiallelic++; missing_haps = num_samples;
			}
			
			// NOW WE HAVE VALID VARIANT,  THAT IS BIALLELIC, SEGREGATING, SNP, WITH NONMISSING ANC AND REF
			else {
				// adjust major/minor allele according to ancestral. // (does nothing unless ancestral is different from ref, then takes complement of non segregating ( we are handling multiallelic case elsewhere))
				if(anc_str.compareToIgnoreCase(snp.alt_allele)==0) {non_ancestral_haps = (num_samples - missing_haps) - non_ancestral_haps;}
			
				// check if ANC synchs with VCF ref/alt as defined by SNP_unit
				else if(anc_str.compareToIgnoreCase(snp.ref_allele)!=0) {
					sites_mismatch_from_anc++; missing_haps = num_samples;
				}
				
				// check full vs partial missing status
				if(snp.full_missing) {snps_missing_from_vcf++;}
				else if(snp.part_missing) { snps_missing_partial_from_vcf++; } 
				else {}
				
			}

			
		}
		
		
		/////////////////////////////////////////////////////

		
		// quick check to make sure SNP position was VALID or halt
		if(variant != null && variant.getStart() < pos) { // the SNP has not been updated and lags behind the position 
			System.out.println("snp needs to be ahead of current position"); System.exit(0);
		}
		
		
		
		// returns 2 numbers. the number of missing sites, and the number of non-ancestral sites.
		return new int[] {missing_haps, non_ancestral_haps, info_at_position};
	}

	
	
	///////////////////////////////////////////////////////////////
	// METALOCUS STREAM CONVERSION //////
	/////////////////////////////////////
	// return reduced stream if mL_size=1 otherwise convert to mL stream
	
	int[][] getStream(int mL_size) throws FileNotFoundException{
	
		// if mL size is 1, just return regular locus skipping reducedStream
		int[][] rS = this.getReducedStream();	
		if(mL_size == 1) {return rS;} // locus skipping, reduced stream for this case. 
		
		// otherwise use reducedStream to create the metaLocus stream

		int num_mL = (this.chrom_seg_rb - this.chrom_seg_lb)/mL_size + 1 ;
		
		int[][] out = new int[num_samples][num_mL];
		int mL_i = 0;
		int mL_remainder = mL_size;
		
		
		int rS_i = 0;
		int[] tract = new int[] {rS[0][rS_i], rS[1][rS_i], rS[2][rS_i]};
		// tract is counting tuple of (head_#ss , tract_remainder, type )
		
		
		while(rS_i < rS[0].length) {
			
			// how far along chromosome to progress 
			// (after processing head of tract)
			int distance = Math.min(mL_remainder, tract[1]);
			
			
			// process head of tract
			// then the distance till new mL or tract
			
			if(tract[2]==1) { // if ss tract, this contributes to mL
				
				//process head
				out[tract[0]][mL_i]++; 
				// alter head to nss in case this tract is hanging later
				tract[0] = 0; 				
				
				// process tail
				out[0][mL_i]=out[0][mL_i] + distance - 1;
				
			}
			
			else if(tract[2] == 2) {} // no contribution for missing data
			
			else {System.err.print("Invalid Tract Type"); System.exit(0);}
			
			// shorten the remainders appropriately.
			// at least one should be zero
			mL_remainder = mL_remainder - distance;
			tract[1] = tract[1] - distance;
		

			// CLOSE AND OPEN next mL or tract if needed
			
			if(mL_remainder == 0) { // close mL
				mL_i++; 
				mL_remainder = mL_size;
			}
					
			if(tract[1] == 0) { // close tract
				rS_i++; 
				
				// make sure havent hit end, then get next tract
				if(rS_i < rS[0].length) { 
					tract = new int[] {rS[0][rS_i], rS[1][rS_i], rS[2][rS_i]};
				}
				
			}
					
			
			
		}
			
		
		return out;
	}
	
	///////////////////////////////////////////////////////////////
	// AUXILIARY METHODS + CLASSES ////////////
	////////////////////////////////////////////
	private boolean isATCG(char c) {
		c = Character.toUpperCase(c);
				
		if(c=='A' || c=='T' || c=='C' || c == 'G') {return true;}
		
		return false;
		
	}
	

	private String get_stats_string() {
		String out = "";
		out = out + "-----------------------------------------------------" + "\n";
		out = out + "**** INFORMATIVE SNPS: " + valid_snps + " ****" + "\n";
		out = out + "Missing Loci in Ref: " + sites_missing_from_ref + "\n";
		out = out + "Total Variant Lines Processed from VCF: " + vcf_variant_count + "\n";
		out = out + "VCF Variants Skipped due to Degenerate Positions: " + vcf_duplicate_position + "\n";
		out = out + "Variants in VCF, missing Ref: " + sites_missing_ref_notVcf + "\n \n";
		
		out = out + "Non-SNPs: " + non_snps + "\n";
		out = out + "SNPs Missing from Anc: " + snps_missing_from_ancestral + "\n";
		out = out + "SNPs Fully Missing in VCF: " + snps_missing_from_vcf + "\n";
		out = out + "SNPS Partially Missing in VCF: " + snps_missing_partial_from_vcf + "\n";
		out = out + "SNPS w/ Multi-Alt Alleles: " + this.snps_multiallelic + "\n";
		out = out + "-----------------------------------------------------" + "\n";
		
		return out;
	}
	
	
	void print_stream_to_file(String out_file) {
		try {
			FileWriter scribe = new FileWriter(out_file);
			int[][] stream = getReducedStream();
			int position = chrom_seg_lb;
		
			scribe.append("PROCESSED " + this.target_vcf_file + "\n");
			scribe.append( "LOCI: " + "FROM " + chrom_seg_lb + " TO " + chrom_seg_rb + "\n");	
			scribe.append("SAMPLE HAPLOTYPE INDICES: " + Arrays.toString(hap_indices) + "\n \n");
			
			scribe.append(get_stats_string() + "\n");	
			scribe.append("tract_position,  number_seg_sites,  tract_type,  tract_length, " + "\n");
			for(int i = 0 ; i < stream[0].length ; i++) {
				
				scribe.append(position + ", " + stream[0][i] + ", " + stream[2][i] + ", " + stream[1][i] + "\n");
				position = position + stream[1][i];

			}
			
			
			scribe.flush();
			scribe.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	
	// supplemental method, static,  obtain length of chromosome from ref file, and check it matches anc file.
	static private int check_lengths(String r, String a){
		
		FastaStream ref = new FastaStream(r,1);
		FastaStream anc = new FastaStream(a,1);
		
		int cl = 0; // track length of ref file in bases. assume only 1 cromosomal segment
		
		char nc_r = ref.next_char();
		char nc_a = anc.next_char();
		
		while(nc_r != (char) 0 && nc_a != (char) 0) {
			nc_r = ref.next_char();
			nc_a = anc.next_char();
			cl++ ;
		}
		
		if(nc_r != nc_a) {
			System.err.println("Lengths of ancestral and reference files are not synched."); System.exit(0);
		}
		
		
		
		return cl;
	}
	
	static private final class FastaStream {

	    private String description;
	    private BufferedReader buff_read;
	    
	    public FastaStream(String filename, int start_site)
	    {	
	    	// setup the buffer, and bring it to the start of the sequence with the description read already
	    	try {
				buff_read = new BufferedReader(new FileReader(filename));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	       	
	       	try {
	       		
	       		int first_char_i  = buff_read.read();
	       		if( first_char_i == -1 ) { throw new IOException( filename + " is an empty file" ); }
	   	     
	       		// now find first char thats not newline, or space
	   
	     		
	    		if( (char) first_char_i != '>' ) {
	    			throw new IOException( "First line of " + filename + " should start with '>'" );
	    		}
	       		description = buff_read.readLine();
	    		
	    		// now we have read first line of file, the description for sequence
	       		
	       	}
	       	catch (IOException e) {
   				// TODO Auto-generated catch block
   				e.printStackTrace();
   			}
	       	
	       	//NOW buffer should be ready after the entire first line. 
	       	
	       	// fastforward to desired start position, so that reading the next base will give the start positions'
	       	if ( start_site < 1 ) {System.out.println("Invalid start position"); System.exit(0);}
	       	char s;
	       	for(int i = 1 ; i < start_site ; i++) {s = next_char();}
	       	
	    }

	   
	    
	    public char next_char() {
	    	int next_char_i = 0;
			try {
				next_char_i = buff_read.read();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	    	
			// if end of stream, return null char
	    	if(next_char_i == -1) {return (char) 0;}
	    	
	    	
	    	char next_char = (char) next_char_i;
	    	
	    	// for newline characters or space characters, ignore and keep going
	    	if(next_char == '\n' || next_char == '\r' || next_char == ' ') {return next_char();}

	    	// end of current sequence, so return null char
	    	if(next_char == '>') {return (char) 0;}
	    	
	    	
	    	return next_char;
	    	
	    }
	    
	    //return first description as String
	    public String getDescription(){
	    	return description;
	    }

	  
	  
	}
	
	
	
	
	private final class SNP_unit{
		// goal is to have this class extract necessary information for a variantContext that we know is a snp.
		// it will process only haplotypes indexed by hap_indices, and pull relevant info from that.
		// in this indexing, 1= individual1_hap1, 2= individual1_hap2, 3=individual2_hap1, 4 = individual2_hap2, ... etc.
		// 
		
		private VariantContext snp;
		String ref_allele;
		String alt_allele;
		
		// all the following stats for the SNP are computed in the relevant sample range only.
		private int missing_haps;
		private int der_haps; 
		
		boolean full_missing = false;	 // fully missing in the relevant sample range
		boolean part_missing = false;	 // partially missing in the relevant sample range
		boolean no_alternate = false;    // non-seg in relevant sample range
		boolean is_multiallelic = false; // marks it being multiallelic in the relevant sample range
		
		SNP_unit(VariantContext v){
			snp = v;
			ref_allele = ""; 
			alt_allele = ""; 
			missing_haps = 0;
			der_haps = 0;

			
			int types_of_alleles = 0;
		
			int num_processed_haps = 0;
			int hap_index = 0;
			GenotypesContext geno_list = snp.getGenotypes();
			
			// double for loop to go through each allele in entire vcf, one individual at a time, and each allele for the individual at a time
			for(int g = 0 ; g < geno_list.size(); g++) {
				List<Allele> allele_list = geno_list.get(g).getAlleles();
				for(int a = 0 ; a < allele_list.size() ; a++) { // iterate through alleles of each individual's genotype
					
					hap_index++; // this reflects the allele index we're on (indexed by haplotype starting with 1)
					
					if(spec_haps_contain(hap_index)) {	// check to make sure hap is one of the specified ones
						Allele allele = allele_list.get(a);
						if(allele.getBaseString().compareTo(".")==0) {missing_haps++;}					
						else {
							// record the first non-missing allele as reference
							if(types_of_alleles == 0) { ref_allele = allele.getBaseString(); types_of_alleles++; alt_allele = allele.getBaseString();}
							
							// record second non-missing type as alternate
							else if(types_of_alleles == 1) { String current_allele = allele.getBaseString();
														
								// check if we consider this one alternate
								if(current_allele.compareToIgnoreCase(ref_allele) != 0 ) { 
									alt_allele = current_allele;  types_of_alleles++;
									der_haps++;
								}
								// otherwise it is same as reference, and we don't count anything
							}
							
							// conk out if we find another type of allele
							else if(types_of_alleles == 2){ String current_allele = allele.getBaseString();
								if(current_allele.compareToIgnoreCase(ref_allele) == 0) {}
								else if(current_allele.compareToIgnoreCase(alt_allele) == 0) {der_haps++; }
								else {	is_multiallelic = true; missing_haps = -1; der_haps = -1; types_of_alleles = 3;}
							}
							// do nothing, already set appropriate markers ^
							else {
								// should never go here
								assert(false);
							}
							
						}
						num_processed_haps++; // counter for number of haps we processed
					}
					// finished processing the single allele for the hap that is in valid interval
				}
			}
			
			//check to make sure that we processed all num_samples haplotypes
			if(num_processed_haps != num_samples) { System.out.print("Specified hap is out of range of VCF, or index used multiple times in one group. Try adjusting the haplotype grouping."); System.exit(0);}
			
			// tag certain properties based on read info
			if (missing_haps > 0 ) {part_missing = true;}
			if (types_of_alleles == 0) {full_missing = true; part_missing = false;}
			if (types_of_alleles == 1) {no_alternate = true;}
			
			
		}
		
		int getMissing() {return missing_haps;}
		
		int getDerived() {return der_haps;}
		
		
		
		
		
	}
	
	private boolean spec_haps_contain(int x) {
		for(int i = 0 ; i < hap_indices.length; i++) {
			if(hap_indices[i] == x) {return true;}
		}
		return false;
	}
		
	void update_max_tract_size(int ms) {
		max_tract_size = ms;
	}
	

}
