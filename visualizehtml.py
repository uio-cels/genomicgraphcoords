from __future__ import absolute_import
from __future__ import division
import numpy as np

DEBUG = False

class VisualizeHtml(object):
    """
    Attempt to make a simple html visualization
    """

    def __init__(self, graph, minOffset, maxOffset, id, levels, description='', width=800, genes=[], start_position=None, trans=None):


        self.padding = 50  # Gap between blocks
        self.gap_pixels = 0  # Extra gap pixels that are added
        self.graph = graph
        self.color_counter = 4
        self.colors = ["#0000aa", "#5F96C5", "#C58E8E", "#cccccc", "#6699ff", "orange", "indigo"]
        self.gene_colors = ["darkorange", "#D9EDF7", "#aaaaaa", "pink", "#6699ff",
                            "orange", "#00ffff", "#99ffcc", "darkgray",
                            "#ffff00", "#999966",
                            "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray",
                            "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray",
                            "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray"
                            ]
        self.gene_counter = 0
        self.genes_plotted_heights = {} # Dict of heights for genes
        self.offset_positions = {}  # Dict of offset pos for each alt loci/chrom
        self.offset_counter = 0
        self.lowest_offset = 0
        self.vis_id = id
        self.genes = genes #list(reversed(intervals))
        self.levels = levels
        self.trans = trans
        self.start_position = start_position

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



        # Produce gene selection html

        if len(genes) > 0:
            html_gene_selection = "<p><b>Check the genes that should be displayed (maximum 3)</b> <span style='margin-left: 20px;' id='gene_selector_message'></span></p>"

            for i, gene in enumerate(self.genes):
                checked_status = ""
                if i <= 2:
                    checked_status = "checked"
                html_gene_selection += """

                <span style='display: inline-block; padding: 3px;'>
                <label style='font-weight: 100;'>
                    <input type='checkbox' id='checkbox_%d' onclick='show_gene(%d);' %s> %s
                </label>
                </span>
                """ % (i, i, checked_status, gene.name)
        else:
            html_gene_selection = "<p>No genes in this area</p>"

        self.html += """
        <div class='row'>
            <div class='col-md-8'>
                <h4>%s</h4>
            </div>
            <div class='col-md-4'>
				<!--<div class='options'>
					<p><label><input type='checkbox' onclick="$('.exon').toggle();"> Click to show exons <span id='exon_cnt'></span></label></p>
				</div>
			    -->
			</div>
        """ % description

        self.html += """
            <div class='col-md-12'>
        """

        self.html += html_gene_selection

        self.html += """
                <div class='visualization'>
                    <div style='position: relative;
                                float: right;
                                width: 150px;
                                background-color: white;
                                height: 20px; margin-top: 5px;'>
                        <p style='font-size: 0.8em;'>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Main path (GRCh38)</font><br>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Flanking regions </font><br>
                            <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span> <font color='black'>Alternative locus</font>
                        </p>
                    </div>
        """ % (self.colors[3], self.colors[2], self.colors[1])

        self.html_blocks = {}  # Dict with blocks as keys. Includes only html for block
        self.html_intervals = {}  # Html for all genes. Key is blocks and interval number
        self.html_exons = {}   # Neste dict, key is block and interval
        for b in self.graph.blocks:
            self.html_intervals[b] = {}
            self.html_exons[b] = {}

            for i, gene in enumerate(genes):
                self.html_exons[b][i] = ""
                self.html_exons[b][i] = ""


        # Gene labels
        self.html += """
        <div style='position: relative;
                    float: left;
                    width: 400px;
                    background-color: white;
                    height: 20px; margin-top: 5px; margin-left: 10px'>

        """
        i = 0

        for gene in genes:

            display = "None"
            if i < 3:
                display = "block"

            self.html += """
                <div id='label_%d' class='interval-label' style='display: %s; font-size: 0.8em;'>
                <span style='background-color: %s; width: 30px; height: 12px; display: inline-block'></span>
                 <font color='black'>%s</font><br>
                 </div>
                """ % (i, display, self.gene_colors[i%len(self.gene_colors)], "Gene: " + gene.name + " (" + gene.name + ")")

            if i >= 3:
                # Print javascript for hiding interval
                self.html += """
                <script>
                $(document).ready(function(){
                    $('.interval_%s').hide();
                });
                </script>
                """ % (i)


            i += 1



        self.html += """

        </div>
        """


        self.visualize_v2()
        self.visualize_genes()
        self.html += self.html_arrows

        # Write js to set number of exons
        if id == 0:
            self.html += """
            <script>
                $(document).ready(function(){
                    $('#exon_cnt').html('(%d)');
                });
            </script>""" % self.exon_cnt



    def visualize_genes(self):
        """
        Visualizes genes
        :param intervals: A list of Genes (of type offsetbasedgraph.graphutils.Gene)
        """
        for gene in self.genes:
            #print("<p>Visualizing gene</p>")
            #print("gene")

            interval = gene.transcription_region
            self._plot_interval(interval, gene.name)
            for exon in gene.exons:
                self._plot_exon_interval(exon, interval)


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


    def _plot_exon_interval(self, interval, parent_interval):
        self.exon_cnt += 1
        block = interval.region_paths[0]
        #print(interval)
        if block not in self.graph.blocks:
            return

        assert(block in self.graph.blocks)
        # Find parent interval start in this
        parent_start = 0
        if parent_interval.start_position.region_path_id == block:
            parent_start = parent_interval.start_position.offset

        start = (interval.start_position.offset - parent_start) * self.width_ratio
        end = (interval.end_position.offset - parent_start) * self.width_ratio

        self._plot_interval_in_block(start, end, 0, interval, block, 0, "Exon", True)

    def _plot_interval(self, interval, name, is_exon = False):

        if is_exon:
            self.exon_cnt += 1
            #print("<p><b>Plotting exon: %s</b></p>" % interval)

        for block in interval.region_paths:
            if not block in self.block_positions:
                if DEBUG: print("Warning. Block %s not found in offset_positions when visualizing interval" % (block))
                continue

            #plot_info = self.offset_positions[block]
            pos = self.block_positions[block]
            start = 0 #pos[0]
            #print("<p>Start pos is %d</p>" % start)
            #end = pos[0] + pos[2]
            end = pos[2]
            if block == interval.region_paths[0]:
                #start += interval.start_position.offset * self.width_ratio
                start = interval.start_position.offset * self.width_ratio
            if block == interval.end_position.region_path_id:
                #end = pos[0] + interval.end_position.offset * self.width_ratio
                end = interval.end_position.offset * self.width_ratio

            parent_width = pos[2]

            if is_exon:
                self._plot_interval_in_block(start, end, pos[1], interval, block, parent_width, name, True)
            else:
                self._plot_interval_in_block(start, end, pos[1], interval, block, parent_width, name)



    def _plot_interval_in_block(self, start, end, level, interval_obj, block, parent_width, name = "", is_exon = False):

        #level += 0.3 + 0.3 * (self.color_counter - 4)
        #print "=
        if end - start == 0:
            return

        top = 1 + self.gene_height * self.gene_counter
        color = self.gene_colors[self.gene_counter]
        height = self.gene_height
        classname = "interval"
        classname2 = "interval_%d" % self.gene_counter
        margin_right = parent_width - end
        position = "relative"

        if is_exon:
            classname = "exon"
            classname2 = "exon_%d" % self.gene_counter
            height = self.exon_height
            color = "black"
            top = 1 #(top + (self.gene_height - self.exon_height) /  2.0)
            margin_right = 0
            position = "absolute"


        html = ""
        html += "<div class='%s %s'" % (classname, classname2)

        html += " style='z-index: 10;"
        #html += "left: %.2fpx;" % start
        html += "left: %.2fpx;" % start
        html += "margin-right: %.2fpx;" % margin_right
        html += "width: %.2fpx;" % (max(1, end - start))
        if is_exon:
            html += "top: %.2fpx;" % (top)

        html += "height: %dpx;" % (height)
        html += "background-color: %s;" % color
        html += "position: %s;" % position
        if not is_exon:
            html += "display: block;"
        else:
            html += "display: inline-block;"

        html += "' "
        html += "data-parent-width='%d'" % parent_width
        if not is_exon:
            html += "data-interval-id='%d'" % self.gene_counter

        html += "data-notation='%s'" % interval_obj.notation()
        html += "data-gene-name='%s'" % name
        html += "data-gene-name2='%s'" % name
        html += "data-graph-id='%d'>" % self.vis_id

        self.genes_plotted_heights[name] = self.gene_counter

        if not is_exon:
            self.html_intervals[block][self.gene_counter] = html
        else:
            #print(self.gene_counter)
            #print(self.html_exons[block])
            self.html_exons[block][self.gene_counter] += html + "</div>"


    def _coordinate(self, rp):
        """
        Returns the hierarhcial and sequential coordinates of a region path
        """

        length = self.graph.blocks[rp].length()
        # Translate rp back to get GRCh38 hier. coordinates
        from offsetbasedgraph import Interval, Graph

        hier_id = str(rp)
        hier_of = 0

        origin = Graph.block_origin(rp)
        if origin == "main" or origin == "merged":
            dist_back = self._distance_to_start(rp)
            hier_id = self.start_position.region_path_id
            hier_of = dist_back + self.start_position.offset

        return (str(rp), "0", str(hier_id), str(hier_of), str(length))

    def _plot(self, xstart, xend, level, color, rp_id):

        html = ""

        y = self.block_height * 2 * ( level + 1)
        x = self.gap_pixels + (xstart - self.minOffset) * self.width_ratio
        width = (xend - xstart) * self.width_ratio

        html += "<div class='block' style='position: absolute;"
        html += "left: %.2fpx;" % x
        html += "width: %.2fpx;" % width
        html += "top: %.2fpx;" % (y)
        html += "height: %dpx;" % (self.block_height)
        html += "background-color: %s;" % color
        html += "' "
        html += " data-rpid='%s'" % (rp_id)
        html += " data-rpname='%s'" % (self._pretty_alt_loci_name(rp_id))
        html += " data-graph-id='%d'" % (self.vis_id)
        html += " data-coordinate='%s'" % ','.join(self._coordinate(rp_id))
        html += ">"

        return x + width, y, width, x, html


    def _plot_level(self, block):
        # Finds the level. Rule: if block contains more than one linea reference
        # it is shared, else, find out whether it is alt or consensus

        return self.levels[block]

        """
        if len(list(block.linear_references.values())) > 1:
            #if DEBUG: print "MOre than one linear reference ::::::::"
            return 1
        else:
            if "alt" in list(block.linear_references.values())[0].chromosome:
                return 0
            else:
                return 2
        """

    def _pretty_alt_loci_name(self, id):
        return id
        return self.graph.pretty_alt_loci_name(id)



    def _scale(self, n):
        return n
        return np.log(n+1)

    def _plot_arrow(self, xstart, ystart, xend, yend):
        """ Plots and arrow
        """
        #print("Plotting arrow  from %d,%d to %d,%d" % (xstart, ystart, xend, yend))
        self.html_arrows += "<div style='position: absolute;"

        if yend < ystart:
            arrow = "short"
            if ystart - yend >= self.block_height * 1:
                arrow = "long"
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart - (ystart - yend) + self.block_height/2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow_up_%s.png' style='" % arrow
            self.html_arrows += "height: %dpx;" % (ystart - yend)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"
        elif yend == ystart:
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart + self.block_height/2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow.png' style='"
            self.html_arrows += "height: %dpx;" % (self.block_height/4)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"
        else:
            arrow = "short"
            if yend - ystart >=  self.block_height * 1:
                arrow = "long"
            if xend - xstart > 100:
                arrow = "wide"
            self.html_arrows += "left: %dpx;" % xstart
            self.html_arrows += "top: %dpx;" % (ystart + self.block_height/2)
            self.html_arrows += "'>"
            self.html_arrows += "<img src='arrow_down_%s.png' style='" % arrow
            self.html_arrows += "height: %dpx;" % (yend-ystart)
            self.html_arrows += "width: %dpx;" % (xend - xstart)
            self.html_arrows += "'>"

        self.html_arrows += "</div>"

    def _get_longest_previous_block(self, block):
        back_block = block
        g = self.graph
        if len(g.reverse_adj_list[back_block]) > 1 and \
                g.blocks[g.reverse_adj_list[back_block][1]].length() > \
                g.blocks[g.reverse_adj_list[back_block][0]].length():
            back_block = g.reverse_adj_list[back_block][1]
        else:
            back_block = g.reverse_adj_list[back_block][0]

        return back_block

    def _distance_to_start(self, b):
        # Find distance back to start
        g = self.graph

        if b == g.start_block:
            return 0
        back_block = self._get_longest_previous_block(b)
        #print("Finding back for %s" % b)
        distance = 0
        while True:
            block_size = self._scale(g.blocks[back_block].length())
            distance += self.padding / self.width_ratio + block_size
            if back_block == g.start_block:
                break

            # Choose longest path back
            back_block = self._get_longest_previous_block(back_block)

        return distance

    def visualize_v2(self):
        # Try to visualize more complex graphs
        block = self.graph.start_block #self.graph.blocks[self.graph.start_block]
        #self.offset_counter = self._scale(list(block.linear_references.values())[0].start)
        self.offset_counter = self._scale(0)


        # Find x position of all blocks
        self.block_positions = {}
        for b in self.graph.blocks:
            start = self._distance_to_start(b)
            end = start + self._scale(self.graph.blocks[b].length())

            # Plot rp
            xend, y, width, xstart, html = self._plot(start, end, self.levels[b], self.colors[self.levels[b] + 1], b)
            self.html_blocks[b] = html

            self.block_positions[b] = (xstart, y, width)
            #print("<p>PLotted %s at %d,%d with width %d, ending at %d</p>" % (b, xstart, y, width, xend))


        # Plot arrows using all edges
        g = self.graph
        for b in g.blocks:
            for edge in self.graph.adj_list[b]:
                length = g.blocks[b].length()
                xstart = self.block_positions[b][0] + self._scale(length) * self.width_ratio
                ystart = self.block_positions[b][1]
                xend = self.block_positions[edge][0]
                yend = self.block_positions[edge][1]
                self._plot_arrow(xstart, ystart, xend, yend)


        return

    def __str__(self):
        html = self.html
        for block, block_html in self.html_blocks.items():
            html += block_html
            for interval, interval_html in self.html_intervals[block].items():
                html += interval_html
                exon_html = self.html_exons[block][interval]
                html += exon_html
                html += "</div>"

            html += "</div>"

        html += "</div></div></div>"
        return html

    def get_wrapped_html(self):
        # Html wrapped with js includes etc
        html = """
        <html>
        <head>
            <script   src="https://code.jquery.com/jquery-2.2.4.min.js"   integrity="sha256-BbhdlvQf/xTY9gja0Dq3HiwQF8LaCRTXxZKRutelT44="   crossorigin="anonymous"></script>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">

            <!--<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css" integrity="sha384-fLW2N01lMqjakBkx3l/M9EahuwpSfeNvV63J5ezn3uZzapT0u7EYsXMjQV+0En5r" crossorigin="anonymous">-->

            <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>
        </head>
        <body>
            <div class="container" style='min-height: 800px;'>
        """
        html += self.__str__()
        html += "<br><br><br><br><br><br><br></div></body></html>"

        return html
