package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.alg.CycleDetector;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Efficient algorithm to obtain the list of best haplotypes given the {@link BaseGraph instance}.
 *
 * minos: GATK 4.5.0.0's class, unchanged except that each search reports to MinosKBest the smallest
 * max-num-haplotypes-in-population that reproduces it (counting only; every decision is the stock one).
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GraphBasedKBestHaplotypeFinder<V extends BaseVertex, E extends BaseEdge> extends KBestHaplotypeFinder<V, E> {

    public final Comparator<KBestHaplotype<V, E>> K_BEST_HAPLOTYPE_COMPARATOR = Comparator.comparingDouble(KBestHaplotype<V, E>::score)
            .reversed()
            .thenComparing(KBestHaplotype<V, E>::getBases, BaseUtils.BASES_COMPARATOR.reversed()); // This is an arbitrary deterministic tie breaker.

    /**
     * Constructs a new best haplotypes finder.
     *
     * @param graph the graph to search.
     * @param sources source vertices for all haplotypes.
     * @param sinks sink vertices for all haplotypes.
     *
     * @throws IllegalArgumentException if <ul>
     *     <li>any of {@code graph}, {@code sources} or {@code sinks} is {@code null} or</li>
     *     <li>any of {@code sources}' or any {@code sinks}' member is not a vertex in {@code graph}.</li>
     * </ul>
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph, final Set<V> sources, final Set<V> sinks) {
        super(sinks, sources, graph);
    }

    /**
     * Constructor for the special case of a single source and sink
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph, final V source, final V sink) {
        this(graph, Collections.singleton(source), Collections.singleton(sink));
    }

    /**
     * Constructor for the default case of all sources and sinks
     */
    public GraphBasedKBestHaplotypeFinder(final BaseGraph<V, E> graph) {
        this(graph, graph.getSources(), graph.getSinks());
    }

    @Override
    public boolean keepCycles() {
        return false;
    }

    /**
     * Implement Dijkstra's algorithm as described in https://en.wikipedia.org/wiki/K_shortest_path_routing
     */
    @Override
    public List<KBestHaplotype<V, E>> findBestHaplotypes(final int maxNumberOfHaplotypes) {
        final List<KBestHaplotype<V, E>> result = new ArrayList<>();
        final PriorityQueue<KBestHaplotype<V, E>> queue = new PriorityQueue<>(K_BEST_HAPLOTYPE_COMPARATOR);
        sources.forEach(source -> queue.add(new KBestHaplotype<>(source, graph)));

        final Map<V, MutableInt> vertexCounts = graph.vertexSet().stream()
                .collect(Collectors.toMap(v -> v, v -> new MutableInt(0)));

        int minosMaxCount = -1;              // minos: the largest vertex count compared with the cap
        boolean minosBound = false;          // minos: the cap refused to extend a vertex
        while (!queue.isEmpty() && result.size() < maxNumberOfHaplotypes) {
            final KBestHaplotype<V, E> pathToExtend = queue.poll();
            final V vertexToExtend = pathToExtend.getLastVertex();
            if (sinks.contains(vertexToExtend)) {
                result.add(pathToExtend);
            } else {
                final int minosCount = vertexCounts.get(vertexToExtend).getAndIncrement();
                minosMaxCount = Math.max(minosMaxCount, minosCount);
                minosBound |= minosCount >= maxNumberOfHaplotypes;
                if (minosCount < maxNumberOfHaplotypes) {
                    final Set<E> outgoingEdges = graph.outgoingEdgesOf(vertexToExtend);
                    int totalOutgoingMultiplicity = 0;
                    for (final BaseEdge edge : outgoingEdges) {
                        totalOutgoingMultiplicity += edge.getMultiplicity();
                    }

                    for (final E edge : outgoingEdges) {
                        queue.add(new KBestHaplotype<>(pathToExtend, edge, totalOutgoingMultiplicity));
                    }
                }
            }
        }
        // minos: a search that ran out of paths with no vertex refused is the same search under any cap of at least
        // max(paths found, largest vertex count + 1): both comparisons with the cap give the same answers
        org.broadinstitute.hellbender.utils.MinosKBest.record(Math.max(result.size(), minosMaxCount + 1),
                minosBound || !queue.isEmpty());
        return result;
    }

}
