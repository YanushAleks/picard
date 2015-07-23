package picard.analysis;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;

public class WgsProcessingData implements Runnable {

	private LinkedBlockingQueue<SamLocusIterator.LocusInfo> infoQueue;
	private AtomicLong basesExcludedByBaseq;
	private AtomicLong basesExcludedByOverlap;
	private AtomicLong basesExcludedByCapping;
	
	private AtomicLongArray HistogramArray;
	private AtomicLongArray baseQHistogramArray;
	
	private AtomicLong countNonNReads;
	
	private int MINIMUM_BASE_QUALITY;
	private int max;
	private final ProgressLogger progress;
	
	private boolean usingStopAfter;
	
	private AtomicBoolean haveNewIndex;
	private AtomicLongArray thCurIndex;
	private LinkedList<ReferenceSequence> listRefSeq;
	private final int threadNumber;
	private AtomicLong minIndex;
	
	public WgsProcessingData(LinkedList<ReferenceSequence> listRefSeq, 
			AtomicLongArray thCurIndex, AtomicBoolean haveNewIndex,
			AtomicLong minIndex, 
			LinkedBlockingQueue<SamLocusIterator.LocusInfo> infoQueue, 
			AtomicLong basesExcludedByBaseq, AtomicLong basesExcludedByOverlap,
			AtomicLong basesExcludedByCapping, AtomicLongArray HistogramArray,
			AtomicLongArray baseQHistogramArray, int MINIMUM_BASE_QUALITY, 
			int max, final ProgressLogger progress, AtomicLong countNonNReads, 
			boolean usingStopAfter, int threadNumber) {
		
		this.listRefSeq = listRefSeq;
		this.thCurIndex = thCurIndex;
		this.haveNewIndex = haveNewIndex;
		this.infoQueue = infoQueue;
		this.basesExcludedByBaseq = basesExcludedByBaseq;
		this.basesExcludedByOverlap = basesExcludedByOverlap;
		this.basesExcludedByCapping = basesExcludedByCapping;
		this.HistogramArray = HistogramArray;
		this.baseQHistogramArray = baseQHistogramArray;
		this.MINIMUM_BASE_QUALITY = MINIMUM_BASE_QUALITY;
		this.max = max;
		this.progress = progress;
		this.countNonNReads = countNonNReads;
		this.usingStopAfter = usingStopAfter;
		this.threadNumber = threadNumber;
		this.minIndex = minIndex;
	}
	
	@Override
	public void run() {
		while(true) {
			SamLocusIterator.LocusInfo info = null;
			try {
				info = infoQueue.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
				System.exit(1);
			}
			
            //exit
            if (info.getPosition() == -1) return;
            
            //data processing
            long buferInfoIndex = info.getSequenceIndex();
            long buferIndex = thCurIndex.getAndSet(this.threadNumber, buferInfoIndex);
            if (buferIndex != buferInfoIndex)
            	this.haveNewIndex.set(true);
            final ReferenceSequence ref = this.listRefSeq.get((int)(buferInfoIndex - minIndex.get()));
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (base == 'N') {
                if (this.usingStopAfter)
                	this.countNonNReads.decrementAndGet();
                continue;
            }
            
            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<String>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            
            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {
                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY) {
                    basesExcludedByBaseq.incrementAndGet();
                    continue;
                }
                if (!readNames.add(recs.getRecord().getReadName())) {
                    basesExcludedByOverlap.incrementAndGet();
                    continue;
                }
                pileupSize++;
                if (pileupSize <= max) {
                    baseQHistogramArray.incrementAndGet(recs.getRecord().getBaseQualities()[recs.getOffset()]);
                }
            }
            
            final int depth = Math.min(readNames.size(), max);
            if (depth < readNames.size())
                basesExcludedByCapping.addAndGet(readNames.size() - max);
            
            HistogramArray.incrementAndGet(depth);
            
            progress.record(info.getSequenceName(), info.getPosition());
        }
    }
}
