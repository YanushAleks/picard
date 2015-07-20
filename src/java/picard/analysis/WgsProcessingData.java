package picard.analysis;

import java.util.HashSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;

public class WgsProcessingData implements Runnable {

	private AtomicBoolean endOfRead;
	
	private ReferenceSequenceFileWalker buferRefWalker;
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
	
	public WgsProcessingData(ReferenceSequenceFileWalker refWalker, 
			LinkedBlockingQueue<SamLocusIterator.LocusInfo> infoQueue, 
			AtomicLong basesExcludedByBaseq, AtomicLong basesExcludedByOverlap,
			AtomicLong basesExcludedByCapping, AtomicLongArray HistogramArray,
			AtomicLongArray baseQHistogramArray, int MINIMUM_BASE_QUALITY, 
			int max, final ProgressLogger progress, AtomicLong countNonNReads, 
			AtomicBoolean endOfRead, boolean usingStopAfter) {
		
		this.buferRefWalker = refWalker;
		this.endOfRead = endOfRead;
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
	}
	
	@Override
	public void run() {
		while (!endOfRead.get() || infoQueue.size() != 0) {
			SamLocusIterator.LocusInfo info = null;
			
			try {
				info = infoQueue.poll(50, TimeUnit.MILLISECONDS);
			} catch (InterruptedException e) {
				e.printStackTrace();
				System.exit(1);
			}
            
            //if queue is empty
            if (info == null) continue;
            
            //data processing
            final ReferenceSequence ref = buferRefWalker.get(info.getSequenceIndex());
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
                    //baseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }
            }
            
            final int depth = Math.min(readNames.size(), max);
            if (depth < readNames.size())
                basesExcludedByCapping.addAndGet(readNames.size() - max);
            //basesExcludedByCapping += readNames.size() - max;
            
            //HistogramArray[depth]++;
            HistogramArray.incrementAndGet(depth);
            
            // Record progress and perhaps stop
            progress.record(info.getSequenceName(), info.getPosition());
        }
    }
}