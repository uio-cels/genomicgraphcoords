import matplotlib.pyplot as plt
import numpy as np
import six
from matplotlib import colors
from config import *
from DbWrapper import DbWrapper

class Visualize:

    def __init__(self, graph):
        self.graph = graph
        self.color_counter = 4
        self.colors = list(six.iteritems(colors.cnames))[6:]
        self.colors = ["#000088", "black", "#880000", "#008800", "purple", "orange", "indigo", "black", "black"]
        self.offset_positions = {}  # Dict of offset pos for each alt loci/chrom
        self.offset_counter = 0
        self.lowest_offset = 0
        plt.plot()
        self.visualize()
        plt.ylim(-10, 10)
        self.db = DbWrapper()


    def visualize_intervals(self, intervals):
        """
        Visualizes intervals
        :param intervals: A list of intervals (of type Segment)
        """
        for interval in intervals:
            for block in interval.block_list:
                plot_info = self.offset_positions[block]
                start = plot_info[0]
                end = plot_info[1]
                if block == interval.start_block:
                    start += interval.start_pos
                if block == interval.end_block:
                    end = start + interval.end_pos

                self._plot_interval_in_block(start, end, plot_info[2])

            self.color_counter += 1

    def _plot_interval_in_block(self, start, end, level):
        level += 0.3 + 0.3 * (self.color_counter - 4)

        plt.plot([start, end], [level, level],
                    self.colors[self.color_counter],
                    linestyle = '-',
                    linewidth=2, solid_capstyle="butt")


    def show(self):
        plt.show()

    def _plot_region_path(self, rp, level=0):
        """
        :return:
        """
        start = self.offset_counter
        lr = rp.linear_references.values()[0]
        length = self._scale(lr.end - lr.start)
        #if DEBUG: print "PLotting " + str(rp.id) + "  on level %.2f, %.2f to %.2f" % (level, start, start + length)
        plt.plot([start, start + length], [level, level],
                 color=self.colors[level + 1], linewidth=10, solid_capstyle="butt")

        self.offset_positions[rp.id] = [start, start + length, level]

        return start + length

    def _plot_level(self, block):
        # Finds the level. Rule: if block contains more than one linea reference
        # it is shared, else, find out whether it is alt or consensus
        if len(block.linear_references.values()) > 1:
            #if DEBUG: print "MOre than one linear reference ::::::::"
            return 0
        else:
            if "alt" in block.linear_references.values()[0].chromosome:
                return 1
            else:
                return -1


    def _scale(self, n):
        return n
        return np.log(n+1)

    def visualize(self):
        """
        Simple function for visualizing the graph using line plots in
        matplotlib
        """
        block = self.graph.blocks[self.graph.start_block]
        self.offset_counter = self._scale(block.linear_references.values()[0].start)
        self.offset_counter = self._plot_region_path(block, -1)

        while True:
            #if DEBUG: print "Iteration: " + str(block.id)
            # Plot the next blocks (will never be more than two blocks)
            offsets = []
            for next in self.graph.block_edges[block.id]:
                level = self._plot_level(self.graph.blocks[next])
                next_block = self.graph.blocks[next]
                o = self._plot_region_path(next_block, level)
                block = next_block
                offsets.append(o)

            if block.id not in self.graph.block_edges:
                return
            if len(self.graph.block_edges[block.id]) == 0:
                return

            self.offset_counter = max(offsets)


class VisualizeHtml():
    """
    Attempt to make a simple html visualization
    """

    def __init__(self, graph, minOffset, maxOffset, id, description='', width=800, intervals=[]):

        self.gap_pixels = 0  # Extra gap pixels that are added
        self.graph = graph
        self.color_counter = 4
        self.colors = list(six.iteritems(colors.cnames))[6:]
        self.colors = ["#0000aa", "#5F96C5", "#C58E8E", "#cccccc", "purple", "orange", "indigo"]
        self.gene_colors = ["darkorange", "#D9EDF7", "#aaaaaa", "pink", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"]
        self.gene_counter = 0
        self.genes_plotted_heights = {} # Dict of heights for genes
        self.offset_positions = {}  # Dict of offset pos for each alt loci/chrom
        self.offset_counter = 0
        self.lowest_offset = 0
        self.vis_id = id
        self.intervals = intervals #list(reversed(intervals))

        self.width = width
        self.maxOffset = maxOffset
        self.minOffset = minOffset
        self.width_ratio = float(self.width) / (self.maxOffset - self.minOffset)

        self.width_used = self.width

        self.gene_height = 10
        self.block_height = 35
        self.exon_height = 6

        self.exon_cnt = 0

        self.svg_lines = ""
        self.html_arrows = ""
        self.html = ""
        self.html += """
        <div class='row'>
            <div class='col-md-8'>
                <h4>%s</h4>
            </div>
            <div class='col-md-4'>            
				<div class='options'>
					<p><label><input type='checkbox' onclick="$('.exon').toggle();"> Click to show exons <span id='exon_cnt'></span></label></p>
				</div>
			</div>
        """ % description

        self.html += """
            <div class='col-md-12'>
                <div class='visualization'>
                    <div style='position: relative;
                                float: right;
                                width: 150px;
                                background-color: white;
                                height: 20px; margin-top: 5px;'>
                        <p style='font-size: 0.8em;'>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Main path (GRCh38)</font><br>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Merged </font><br>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Alternative loci</font>
                        </p>
                    </div>
        """ % (self.colors[3], self.colors[2], self.colors[1])


        # Gene labels
        self.html += """
        <div style='position: relative;
                    float: left;
                    width: 400px;
                    background-color: white;
                    height: 20px; margin-top: 5px; margin-left: 10px'>
            <p style='font-size: 0.8em;'>
        """
        i = 0
        for gene in intervals:
            if not gene.is_exon:
                self.html += """
                <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span>
                 <font color='black'>%s</font><br>
                """ % (self.gene_colors[i], "Gene: " + gene.name + " (" + gene.gene_name + ")")
            i += 1


        self.html += """
            </p>
        </div>
        """


        self.visualize()
        self.visualize_intervals()
        self.html += self.html_arrows

        # Write js to set number of exons
        if id == 0:
            self.html += """
            <script>
                $(document).ready(function(){
                    $('#exon_cnt').html('(%d)');
                });
            </script>""" % self.exon_cnt

        self.html += "</div></div></div>"

    def visualize_intervals(self):
        """
        Visualizes intervals
        :param intervals: A list of intervals (of type Segment)
        """
        for interval in self.intervals:
            for block in interval.block_list:
                plot_info = self.offset_positions[block]
                start = plot_info[0]
                end = plot_info[0] + plot_info[2]
                if block == interval.start_block:
                    start += interval.start_pos * self.width_ratio
                if block == interval.end_block:
                    end = plot_info[0] + interval.end_pos * self.width_ratio

                if interval.is_exon:
                    self._plot_exon(start, end, plot_info[1], interval)
                else:
                    self._plot_interval_in_block(start, end, plot_info[1], interval)

            self.color_counter += 1
            self.gene_counter += 1

    def _interval_js_css(self):
        self.html += """
        <style type='text/css'>
            .interval_%d:hover .interval_%d{
                border: 1px solid white;
            }
        </style>
        """ % (self.gene_counter, self.gene_counter)

    def _plot_exon(self, start, end, level, interval_obj):
        image = "exon_start"
        top = level + 1 + self.gene_height * self.genes_plotted_heights[interval_obj.gene_name]

        self.html += "<div class='exon exon_%d'" % self.gene_counter
        self.html += " style='z-index: 12; position: absolute;"
        self.html += "left: %.2fpx;" % (start)
        self.html += "width: %.2fpx;" % (max(2, end - start))
        self.html += "top: %.2fpx;" % (top + (self.gene_height - self.exon_height) / 2.0)
        self.html += "height: %dpx;" % (self.exon_height)
        self.html += "background-color: black;"
        self.html += "' "
        self.html += "data-interval-id='%d'" % self.gene_counter
        self.html += "data-notation='%s'" % str(interval_obj)
        self.html += "data-gene-name='%s'" % interval_obj.name
        self.html += "data-gene-name2='%s'" % interval_obj.gene_name
        self.html += "data-graph-id='%d'>" % self.vis_id
        #self.html += "<img src='%s.png'>" % image
        self.html += "</div>"

        if start < 900:
            self.exon_cnt += 1


    def _plot_interval_in_block(self, start, end, level, interval_obj):
        #level += 0.3 + 0.3 * (self.color_counter - 4)
        #print "=
        if end - start == 0:
            return

        top = level + 1 + self.gene_height * self.gene_counter

        self.html += "<div class='interval interval_%d'" % self.gene_counter

        self.html += " style='z-index: 10; position: absolute;"
        self.html += "left: %.2fpx;" % start
        self.html += "width: %.2fpx;" % (end - start)
        self.html += "top: %.2fpx;" % (top)
        self.html += "height: %dpx;" % (self.gene_height)
        self.html += "background-color: %s;" % self.gene_colors[self.gene_counter]
        self.html += "' "
        self.html += "data-interval-id='%d'" % self.gene_counter
        self.html += "data-notation='%s'" % str(interval_obj)
        self.html += "data-gene-name='%s'" % interval_obj.name
        self.html += "data-gene-name2='%s'" % interval_obj.gene_name
        self.html += "data-graph-id='%d'></div>" % self.vis_id

        self.genes_plotted_heights[interval_obj.name] = self.gene_counter

        #print "%d, %d, %d" % (start, end, level)



        #plt.plot([start, end], [level, level],
        #            self.colors[self.color_counter],
        #            linestyle = '-',
        #            linewidth=2, solid_capstyle="butt")

    def show(self):
        plt.show()

    def _coordinate(self, rp):
        """
        Returns the hierarhcial and sequential coordinates of a region path
        """

        # Sequential coordinates are always id and the first offset is 0
        seqID = rp.id
        seqOf = 0
        # Hierarchical coordinates is same as sequential if this is alternative
        # block. Else, it is the same as on hg38
        if len(rp.linear_references) == 1 and \
                not "alt" in rp.linear_references.values()[0].chromosome:
            hierID = rp.linear_references.values()[0].chromosome
            hierOf = rp.linear_references.values()[0].start
        elif len(rp.linear_references) == 2 and \
            "alt" in rp.linear_references.values()[0].chromosome:
            hierID = rp.linear_references.values()[1].chromosome
            hierOf = rp.linear_references.values()[1].start
        elif len(rp.linear_references) == 2 and \
            "alt" in rp.linear_references.values()[1].chromosome:
            hierID = rp.linear_references.values()[0].chromosome
            hierOf = rp.linear_references.values()[0].start
        else:
            hierID = seqID
            hierOf = seqOf

        # Also get size of region path
        size = max(l.end - l.start for l in rp.linear_references.values())

        return (self._pretty_alt_loci_name(seqID), str(seqOf),\
                self._pretty_alt_loci_name(hierID), str(hierOf), str(size))


    def _plot(self, xstart, xend, level, color, rp):

        y = self.block_height * 2 * ( level + 1)
        x = self.gap_pixels + (xstart - self.minOffset) * self.width_ratio
        width = (xend - xstart) * self.width_ratio

        self.html += "<div class='block' style='position: absolute;"
        self.html += "left: %.2fpx;" % x
        self.html += "width: %.2fpx;" % width
        self.html += "top: %.2fpx;" % (y)
        self.html += "height: %dpx;" % (self.block_height)
        self.html += "background-color: %s;" % color
        self.html += "' "
        self.html += " data-rpid='%s'" % (rp.id)
        self.html += " data-rpname='%s'" % (self._pretty_alt_loci_name(rp.id))
        self.html += " data-graph-id='%d'" % (self.vis_id)
        self.html += " data-coordinate='%s'" % ','.join(self._coordinate(rp))
        self.html += ">"
        self.html += "<font color='white'>%s</font></div>" % ""

        return x + width, y, width, x

    def _plot_region_path(self, rp, level=0):
        """
        :return:
        """
        start = self.offset_counter
        lr = rp.linear_references.values()[0]
        length = self._scale(lr.end - lr.start)
        #if DEBUG: print "PLotting " + str(rp.id) + "  on level %.2f, %.2f to %.2f" % (level, start, start + length)
        xend, y, width, xstart = self._plot(start, start + length, level , self.colors[level + 1], rp)

        self.offset_positions[rp.id] = [xstart, y, width]

        return start + length, xend , y, width

    def _plot_level(self, block):
        # Finds the level. Rule: if block contains more than one linea reference
        # it is shared, else, find out whether it is alt or consensus
        if len(block.linear_references.values()) > 1:
            #if DEBUG: print "MOre than one linear reference ::::::::"
            return 1
        else:
            if "alt" in block.linear_references.values()[0].chromosome:
                return 0
            else:
                return 2

    def _pretty_alt_loci_name(self, id):
        return self.graph.pretty_alt_loci_name(id)



    def _scale(self, n):
        return n
        return np.log(n+1)

    def _plot_arrow(self, xstart, ystart, xend, yend):
        """ Plots and arrow
        """
        #print "Plot from %d,%d to %d,%d" % (xstart, ystart, xend, yend)
        self.html_arrows += "<div style='position: absolute;"

        if yend < ystart:
            arrow = "short"
            if ystart - yend >= self.block_height * 4:
                arrow = "long"
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart - (ystart - yend) + self.block_height / 2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow_up_%s.png' style='" % arrow
            self.html_arrows += "height: %dpx;" % (ystart - yend)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"
        elif yend == ystart:
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart + self.block_height / 2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow.png' style='"
            self.html_arrows += "height: %dpx;" % (self.block_height / 4)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"
        else:
            arrow = "short"
            if ystart - yend >=  + self.block_height * 4:
                arrow = "long"
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart + self.block_height / 2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow_down_%s.png' style='" % arrow
            self.html_arrows += "height: %dpx;" % (yend-ystart)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"

        self.html_arrows += "</div>"


    def visualize(self):
        block = self.graph.blocks[self.graph.start_block]
        self.offset_counter = self._scale(block.linear_references.values()[0].start)
        self.offset_counter, x_coord, y_coord, width = self._plot_region_path(block,
                                                    self._plot_level(block))

        prevEnds = [(x_coord, y_coord, width)]

        while True:
            self.gap_pixels += 20
            #if DEBUG: print "Iteration: " + str(block.id)
            # Plot the next blocks (will never be more than two blocks)
            offsets = []

            ends = []

            for next in self.graph.block_edges[block.id]:
                level = self._plot_level(self.graph.blocks[next])
                next_block = self.graph.blocks[next]
                o, x, y, width = self._plot_region_path(next_block, level)
                block = next_block
                offsets.append(o)
                ends.append((x, y, width))
                self.width_used = max(self.width_used, y + width)

            # Hack
            if width < 1:
                continue

            # PLot arrows from previous end positions (there are 1 or 2),
            # to these end positions (there are 1 or 2)
            for end1 in prevEnds:
                for end2 in ends:
                    self._plot_arrow(end1[0], end1[1], \
                        end2[0]-end2[2], end2[1])

            prevEnds = ends

            if block.id not in self.graph.block_edges:
                return
            if len(self.graph.block_edges[block.id]) == 0:
                return




            self.offset_counter = max(offsets)


    def __str__(self):
        return self.html
if __name__ == "__main__":

    v = Visualize(None)
