package com.actelion.research.chem.hyperspace;

/**
 * Idorsia Pharmaceuticals Ltd. 2020
 * Thomas Liphardt
 *
 * Hyperspace
 */



import com.actelion.research.chem.StructureSearchDataSource;

import java.io.*;
import java.util.*;


/**
 *
 * Class for fast filtering of bitstring fingerprints.
 *
 * Usage:
 *
 * StructureSearchDataSource ssds;
 * BitSetTree bst = BitSetTree.createFromStructureSearchDataSource(ssds, "FragFp" , 512 );
 *
 * to filter for a given FragFp fingerprint and find all rows that pass the filtering:
 *
 * long[] fp;
 * int maxResults = 2000;
 * int rows[] = bst.filterRows(fp, maxResults);
 *
 *
 *
 */
public class BitSetTree implements Serializable {

    private static final long serialVersionUID = 2219385776457978215L;

    //public enum STORAGE_MODE { MEMORY , FILE , ZIP };

    Node root;

    public static class BitSetWithRow {
        public final BitSet bitset;
        public final int row;
        public BitSetWithRow(BitSet bs, int row) {
            this.bitset = bs;
            this.row   = row;
        }
    }




    BitSetTree(Node root) {
        this.root = root;
    }

    /**
     * NOTE!
     * If the recursion finds the first superset in a leaf node,
     * then not all entries of the first_superset node will contain
     * hits! I.e. it is necessary to call collectHits after this!
     *
     *
     * @param b
     * @param first_superset
     * @return

     */
    public boolean testSubset(BitSet b, Node[] first_superset) {
        return this.root.testSubset(b,first_superset);
    }


    public static final class Node implements Serializable {
        // if -1, then this is a leaf
        public int bit;

        BitSet bits_1;
        BitSet bits_0;

        Node left = null;
        Node right = null;
        private List<BitSetWithRow> leaf_data = null;


        public Node(int bit, BitSet bits_0, BitSet bits_1, Node left, Node right, List<BitSetWithRow> data) {
            this.bit = bit;
            this.bits_0 = bits_0;
            this.bits_1 = bits_1;
            this.left = left;
            this.right = right;
            this.leaf_data = data;
        }


        public boolean isLeaf() {
            return this.bit < 0;
        }

        public void setLeft(Node l) {
            this.left = l;
        }

        public void setRight(Node r) {
            this.right = r;
        }

        public int depth() {
            return this.bits_0.cardinality() + this.bits_1.cardinality();
        }


        /**
         * returns the leaf data.
         *
         * @return
         */
        public List<BitSetWithRow> getLeafData() {
            if (this.leaf_data != null) {
                return this.leaf_data;
            } else {

            }
            return null;
        }

        /**
         * NOTE!
         * If the recursion finds the first superset in a leaf node,
         * then not all entries of the first_superset node will contain
         * hits! I.e. it is necessary to call collectHits after this!
         *
         * @param b
         * @param first_superset
         * @return
         */
        public boolean testSubset(BitSet b, Node[] first_superset) {

            first_superset[0] = null;

            if (this.bit < 0) {
                List<BitSetWithRow> n_leaf_data = getLeafData();

                //if(this.leaf_data.isEmpty()) {
                if (n_leaf_data.isEmpty()) {
                    return false;
                }
                boolean subset = true;
                //for(BitSet bsi : leaf_data) {
                for (BitSetWithRow bsi : n_leaf_data) {
                    BitSet ti = (BitSet) bsi.bitset.clone();
                    ti.or(b);
                    if (ti.equals(bsi.bitset)) {
                        first_superset[0] = this;
                        return true;
                    }
                }
                return false;
            }

            BitSet test = (BitSet) bits_1.clone();
            test.or(b);
            if (test.equals(bits_1)) {
                first_superset[0] = this;
                return true;
            }


            if (b.get(this.bit)) {
                return this.right.testSubset(b, first_superset);
            } else {
                return this.right.testSubset(b, first_superset) || this.left.testSubset(b, first_superset);
            }
        }

        /**
         * Returns all rows that contain supersets of the supplied fingerprint.
         *
         * @param fingerprint
         * @param max_results
         * @return
         */
        public int[] filterRows(long[] fingerprint, int max_results) {
            List<BitSetWithRow> results = new ArrayList<>();
            this.collectSuperSetsWithRows( BitSet.valueOf(fingerprint) , results , max_results );
            int[] rows = results.stream().mapToInt(ri -> ri.row).toArray();
            return rows;
        }

        /**
         * Collect supersets of the query bs, including the rows.
         *
         * @param bs
         * @param supersets out-parameter, results will be added to this list
         * @param max_supersets if number of values in the supersets list exceeds this value, the search stops.
         *
         * @return
         */
        public void collectSuperSetsWithRows(BitSet bs, List<BitSetWithRow> supersets, int max_supersets) {
            if(supersets.size() >= max_supersets) {
                return;
            }

            if (this.isLeaf()) { // i.e. we are in a leaf
                //System.out.println("Scan:"+this.getLeafData().size());
                long ts_a1 = System.currentTimeMillis();
                //LinkedList<BitSet> supersets = new LinkedList<>();
                List<BitSetWithRow> n_leaf_data = this.getLeafData();
                //for(BitSet bsi : this.leaf_data) {
                for (BitSetWithRow lei : n_leaf_data) {
                    BitSet ti = (BitSet) lei.bitset.clone();
                    ti.or(bs);
                    if (ti.equals(lei.bitset)) {
                        supersets.add(lei);
                    }
                }
                long ts_a2 = System.currentTimeMillis();
                //System.out.println("T:"+(ts_a2-ts_a1));
            } else {
                if (bs.get(this.bit)) {
                    this.right.collectSuperSetsWithRows(bs, supersets, max_supersets);
                } else {
                    this.left.collectSuperSetsWithRows(bs, supersets, max_supersets);
                    this.right.collectSuperSetsWithRows(bs, supersets, max_supersets);
                }
            }
        }

        /**
         * Collect supersets of the query bs (only the bitsets, without the row).
         *
         * @param bs
         * @return
         */
        public void collectSuperSets(BitSet bs, List<BitSet> supersets) {
            if (this.bit < 0) {
                //LinkedList<BitSet> supersets = new LinkedList<>();
                List<BitSetWithRow> n_leaf_data = this.getLeafData();
                //for(BitSet bsi : this.leaf_data) {
                for (BitSetWithRow lei : n_leaf_data) {
                    BitSet ti = (BitSet) lei.bitset.clone();
                    ti.or(bs);
                    if (ti.equals(lei.bitset)) {
                        supersets.add(lei.bitset);
                    }
                }
            } else {
                if (bs.get(this.bit)) {
                    this.right.collectSuperSets(bs, supersets);
                } else {
                    this.left.collectSuperSets(bs, supersets);
                    this.right.collectSuperSets(bs, supersets);
                }
            }
        }

        public boolean checkAllAreSuperset(BitSet q) {
            if (this.bit < 0) {
                BitSet ti = (BitSet) this.bits_1.clone();
                ti.or(q);
                boolean is_superset = ti.equals(this.bits_1);
                if (is_superset) {
                    return true;
                } // this is only sufficient condition, not necessary.
                // leaf.. check all
                is_superset = true;
                for (BitSetWithRow lei : this.leaf_data) {
                    BitSet bi = lei.bitset;
                    ti = (BitSet) bi.clone();
                    ti.or(q);
                    if (!ti.equals(bi)) {
                        is_superset = false;
                        break;
                    }
                }
                return is_superset;
            }
            if (q.get(this.bit)) {
                if (this.left == null) {
                    return true;
                }
                BitSet empty = new BitSet(q.size());
                // we have to check for any elements in this subtree, therefore the empty bitset.
                if (this.left.getSuperSetIterator(empty).hasNext()) {
                    return false;
                }
                return true;
            }
            return this.right.checkAllAreSuperset(q) && this.left.checkAllAreSuperset(q);
        }


        public SuperSetIterator getSuperSetIterator(BitSet q) {
            return new SuperSetIterator(this, q);
        }

        public int countAll() {
            if (this.isLeaf()) {
                List<BitSetWithRow> n_leaf_data = this.getLeafData();
                //return this.leaf_data.size();
                return n_leaf_data.size();
            }
            int ca = (this.left != null) ? this.left.countAll() : 0;
            int cb = (this.right != null) ? this.right.countAll() : 0;
            return ca + cb;
        }

    }


    /**
     *
     *
     * @param bitsets
     * @param num_bits
     * @param binsize
     * @param max_tries, max. number of bits that are considered to find a balanced bit for splitting the tree into
     *                   two equal parts.
     * @return
     */
    public static BitSetTree createTree( Collection<BitSetWithRow> bitsets , int num_bits , int binsize, int max_tries) {
        Node root = split_recursively(bitsets, new BitSet(num_bits) , new BitSet(num_bits) , num_bits,binsize, "r",max_tries);
        return new BitSetTree(root);
    }

    public static BitSetTree createTree( Collection<BitSetWithRow> bitsets , int num_bits , int binsize) {
        Node root = split_recursively(bitsets, new BitSet(num_bits) , new BitSet(num_bits) , num_bits,binsize, "r");
        return new BitSetTree(root);
    }


    public static Node split_recursively(Collection<BitSetWithRow> bi,BitSet bits_0, BitSet bits_1, int num_bits, int binsize, String tree_pos) {
        return split_recursively(bi,bits_0,bits_1,num_bits,binsize,tree_pos,20);
    }

    public static Node split_recursively(Collection<BitSetWithRow> bi,BitSet bits_0, BitSet bits_1, int num_bits, int binsize, String tree_pos , int max_tries ) {
        if(bi.size() <= binsize) {
            return new Node(-1, bits_0, bits_1, null, null, new ArrayList<>(bi));
        }

        List<Integer> possible_splits = new ArrayList<>();
        for(int zi=0;zi<num_bits;zi++) { if( (!bits_0.get(zi)) && (!bits_1.get(zi)) ) { possible_splits.add(zi); } }
        Collections.shuffle( possible_splits , random );

        double best_split     = 0.0;
        int    best_split_bit = -1;
        for(int zi=0;zi< Math.min( possible_splits.size() , max_tries ) ;zi++)
        {
            //int split = random.nextInt(num_bits);
            int split = possible_splits.get(zi);

            int sa = 0;
            for(BitSetWithRow bsi : bi) {
                sa += (bsi.bitset.get(split))?1:0;
            }
            double split_value = 1.0*sa / bi.size();
            double split_score = Math.min( split_value , 1.0-split_value );
            if( split_score > best_split ) {
                best_split_bit = split;
                best_split     = split_score;
            }
            if( best_split > 0.42 ) {
                break;
            }
        }
        //System.out.println("Node size= " + bi.size() + " Split at bit "+best_split_bit+" , score = " + (best_split) );

        if(best_split_bit<0) {
            System.out.println("wut?");
        }
        else {
            if(true) {
                System.out.println("SplitScore: "+best_split+"  (Size="+bi.size()+",level="+(bits_0.cardinality()+bits_1.cardinality())+")");
            }
        }

        List<BitSetWithRow> bs_a = new ArrayList<>();
        List<BitSetWithRow> bs_b = new ArrayList<>();
        for(BitSetWithRow bsi : bi) {
            if(bsi.bitset.get(best_split_bit)) {
                bs_b.add(bsi);
            }
            else {
                bs_a.add(bsi);
            }
        }

        BitSet bits_0_left  = (BitSet) bits_0.clone();
        BitSet bits_1_left  = (BitSet) bits_1.clone();
        BitSet bits_0_right = (BitSet) bits_0.clone();
        BitSet bits_1_right = (BitSet) bits_1.clone();

        bits_0_left.set(best_split_bit,true);
        bits_1_right.set(best_split_bit,true);

        Node left  = split_recursively(bs_a,bits_0_left,bits_1_left,num_bits,binsize,tree_pos+"_0");
        Node right = split_recursively(bs_b,bits_0_right,bits_1_right,num_bits,binsize,tree_pos+"_1");

        return new Node(best_split_bit,bits_0,bits_1,left,right,null);
    }

    private static Random random = new Random();


    public static boolean isSet(byte[] arr, int bit) {
        int index = bit / 8;  // Get the index of the array for the byte with this bit
        int bitPosition = bit % 8;  // Position of this bit in a byte
        if(index >= arr.length) {
            return false;
        }
        return (arr[index] >> bitPosition & 1) == 1;
    }


    private void collectNodes_dfs(Node n, List<Node> nodes) {
        nodes.add(n);
        if(n.left!=null) {
            collectNodes_dfs(n.left, nodes);
        }
        if(n.right!=null) {
            collectNodes_dfs(n.right, nodes);
        }
    }




    /**
     *
     * @param ssds
     * @param descriptorShortName name of a descriptor that returns long[] objects
     * @param descriptorBits number of bits in the long[] descriptors (i.e. array length * 64)
     * @return
     */
    public static BitSetTree createFromStructureSearchDataSource(StructureSearchDataSource ssds, String descriptorShortName, int descriptorBits) {
        return createFromStructureSearchDataSource(ssds,descriptorShortName,descriptorBits,512,40);
    }

    /**
     *
     * @param ssds
     * @param descriptorShortName
     * @param descriptorBits
     * @param treeBinSize
     * @param maxTries max. number of bits that are considered in each split to find a
     *                 balanced bit for splitting the tree into two equal parts.
     * @return
     */
    public static BitSetTree createFromStructureSearchDataSource(StructureSearchDataSource ssds, String descriptorShortName, int descriptorBits, int treeBinSize, int maxTries) {
        int dcol = ssds.getDescriptorColumn(descriptorShortName);
        List<BitSetWithRow> rows = new ArrayList<>();
        for(int zi=0;zi<ssds.getRowCount();zi++) {
            rows.add(new BitSetWithRow( BitSet.valueOf( (long[]) ssds.getDescriptor(dcol,zi,0, false) ) ,zi));
        }
        BitSetTree bst = BitSetTree.createTree(rows,descriptorBits,treeBinSize,maxTries);
        return bst;
    }


    public static class SuperSetIterator implements Iterator<BitSetWithRow> {

        Node       n                = null;
        BitSet     q                = null;

        Iterator<BitSetWithRow>   current_iterator = null;
        BitSetWithRow                current_next  = null;

        // contains the not yet exhausted nodes (fill in reverse order that you want them to be processed!)
        List<Node>       remaining_childs = null;

        public SuperSetIterator(Node n, BitSet q) {
            this.n    = n;
            this.q    = q;

            // fill remaining_childs (we go "1" first, then "0", therefore add in reverse!):
            if(!n.isLeaf()) {
                this.remaining_childs = new ArrayList<>();
                if(n.left  != null) {this.remaining_childs.add(n.left);}
                if(n.right != null) {this.remaining_childs.add(n.right);}
            }
        }

        @Override
        public boolean hasNext() {
            if(this.current_next==null) {
                tryToFindNext();
            }
            return this.current_next!=null;
        }

        private void tryToFindNext() {
            if(this.current_next!=null) {
                // we already have the next..
                return;
            }

            if(this.n.isLeaf()) {
                if(this.current_iterator == null) {
                    this.current_iterator = this.n.leaf_data.iterator();
                }
                // compute superset tests and check for next:
                while( this.current_iterator.hasNext() ) {
                    //this.current_next = this.current_iterator.next();
                    BitSetWithRow candidate = this.current_iterator.next();
                    // superset test:
                    BitSet ti = (BitSet) candidate.bitset.clone();
                    ti.or(this.q);
                    if(ti.equals(candidate.bitset)) {
                        // ok, we have a next and return
                        this.current_next = candidate;//ti;
                        return;
                    }
                }
                // this leaf is exhausted if we end up here
                this.current_next = null;
                return;
            }
            else {
                // 1. check if we have iterator. if we have, check if we it has next, else set to null
                // 2. if current_iterator == null, check if we have next
                if(this.current_iterator!=null) {
                    if(this.current_iterator.hasNext()) {
                        this.current_next = this.current_iterator.next();
                        return;
                    }
                    else {
                        this.current_iterator = null;
                    }
                }

                while(this.current_next==null) { // && this.remaining_childs.size()>0) {
                    if(this.current_iterator!=null) {
                        if(this.current_iterator.hasNext()) {
                            this.current_next = this.current_iterator.next();
                            return;
                        }
                        else {
                            this.current_iterator = null;
                        }
                    }
                    // we need new iterator?
                    if (this.current_iterator == null) {
                        if (this.remaining_childs.size() > 0) {
                            this.current_iterator = this.remaining_childs.remove(this.remaining_childs.size() - 1).getSuperSetIterator(this.q);
                        } else {
                            this.current_next = null;
                            return;
                        }
                    }
                    // if we end up here, we have next iterator to search through:
                    if(this.current_iterator.hasNext()) {
                        this.current_next = this.current_iterator.next();
                        return;
                    }
                    else {
                        // we continue with this loop..
                        this.current_next = null;
                    }
                }
            }
        }

        @Override
        public BitSetWithRow next() {
            BitSetWithRow next = this.current_next;
            // Set next to null!!
            this.current_next = null;
            return next;
        }
    }




//    public static void main(String args[]) {
//        //test_bitsets();
//        //test_serialization();
//        //test_fetch_supersets();
//        //test_bitsets2();
//        //test_bitsets2_zip();
//        //test_benchmark_a()
//    //    test_benchmark_b();
//    }

    public static void test_bitsets2() {
        Random r = new Random();
        List<BitSetWithRow> test = createRandomBitSetWithRow(r, 20000,64,0.6);

        BitSetTree t  = BitSetTree.createTree( new HashSet<>(test) , 64,8 );

        BitSet to_find = new BitSet(64);
        to_find.set(13);
        to_find.set(17);
        to_find.set(19);
        to_find.set(27);

        List<BitSetWithRow> supersets = new ArrayList<>();
        for(BitSetWithRow ti : test) {
            BitSet tit = (BitSet) ti.bitset.clone();
            tit.or(to_find);
            if(tit.equals(ti.bitset)) {
                supersets.add(ti);
            }
        }

        Node supersets_2[] = new Node[1];
        boolean found = t.testSubset(to_find,supersets_2);

        List<BitSet> supersets_3 = new ArrayList<>();
        t.root.collectSuperSets(to_find, supersets_3);

        if(supersets_3.size()==supersets.size()) {
            System.out.println("ok");
        }
        else {
            System.out.println("not ok");
        }

    }

    public static void test_fetch_supersets() {
        Random r = new Random();
        List<BitSetWithRow> test = createRandomBitSetWithRow(r,64,16,0.6);

        BitSetTree t = BitSetTree.createTree( new HashSet<>(test) , 16,2 );

        BitSet to_find = test.get(3).bitset;
        to_find.set( to_find.nextSetBit(0) , 0 );
        Node supersets[] = new Node[1];
        boolean found = t.testSubset(to_find,supersets);

        //List<BitSet> bs = supersets[0].collectSuperSets(to_find);
        if(found) {
            System.out.println("ok");
        }
        else {
            System.out.println("ERROR.. not working :(");
        }
    }


    public static void test_benchmark_a() {
        Random r = new Random();

        int bits = 2048;
        List<BitSetWithRow> rand_bs = createRandomBitSetWithRow(r,200000,bits, 0.6);

        System.out.println("Create tree:");
        BitSetTree bst = createTree(new HashSet<>(rand_bs),bits,64);
        System.out.println("Create tree: done!");

        //List<BitSetWithRow> rand_test = createRandomBitSetWithRow(r,1000,bits,0.25);
        List<BitSetWithRow> rand_test = createRandomBitSetWithRow(r,1000,bits,0.12);

        System.out.println("Start eval:");
        long ts_a = System.currentTimeMillis();

        for(BitSetWithRow bsi : rand_test) {
            Node[] result = new Node[1];
            System.out.println(bst.testSubset(bsi.bitset,result));
        }

        long ts_b = System.currentTimeMillis();
        System.out.println("Time= "+(ts_b-ts_a));
    }

    public static void test_benchmark_b() {
        Random r = new Random();

        int bits = 512;
        //List<BitSetWithRow> rand_bs = createRandomBitSetWithRow(r,8000000,bits, 0.6);
        List<BitSetWithRow> rand_bs = createRandomBitSetWithRow(r,800000,bits, 0.8);

        System.out.println("Create tree:");
        BitSetTree bst = createTree(new HashSet<>(rand_bs),bits,100000);
        System.out.println("Create tree: done!");

        //List<BitSetWithRow> rand_test = createRandomBitSetWithRow(r,1000,bits,0.25);
        List<BitSetWithRow> rand_test = createRandomBitSetWithRow(r,100,bits,0.3);

        System.out.println("Start eval:");
        long ts_a = System.currentTimeMillis();

        for(int zr=0;zr<rand_test.size();zr++) {
            BitSetWithRow bsi = rand_test.get(zr);
            //Node[] result = new Node[1];
            //System.out.println(bst.testSubset(bsi.bitset,result));

            if(true) {
                long ts_a1 = System.currentTimeMillis();
                SuperSetIterator ssiterator = new SuperSetIterator(bst.root, bsi.bitset);
                List<BitSetWithRow> rows = new ArrayList<>();
                for (int zi = 0; zi < 2000; zi++) {
                    if (ssiterator.hasNext()) {
                        rows.add(ssiterator.next());
                    } else {
                        //System.out.println( zr+" -> Results: "+rows.size());
                        break;
                    }
                }
                long ts_b1 = System.currentTimeMillis();
                System.out.println(zr + " -> Results: " + rows.size());
                System.out.println("Time A: "+ (ts_b1-ts_a1) );
            }
            if(true) {
                long ts_a1 = System.currentTimeMillis();
                List<BitSetWithRow> rows = new ArrayList<>();
                bst.root.collectSuperSetsWithRows(bsi.bitset,rows,2000);
                long ts_b1 = System.currentTimeMillis();
                System.out.println(zr + " -> Results: " + rows.size());
                System.out.println("Time B: "+ (ts_b1-ts_a1) );
            }
        }

        long ts_b = System.currentTimeMillis();
        System.out.println("Time= "+(ts_b-ts_a));
    }



    public static List<BitSetWithRow> createRandomBitSetWithRow(Random r, int n, int bits, double density) {
        List<BitSetWithRow> rand_bs = new ArrayList<>();
        for(int zi=0;zi<n;zi++) {
            BitSet bsi = new BitSet(bits);
            for(int zj=0;zj<bits;zj++) {
                bsi.set( zj , r.nextFloat() < density );
            }
            rand_bs.add( new BitSetWithRow(bsi,zi) );
        }
        return rand_bs;
    }

}
