/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Plotting of data in PDF format
  \author Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifdef HAVE_HARU

#include <itpp/graphics/plot_pdf.h>
#include <itpp/base/itfile.h>


namespace itpp {

  PlotPDF::PlotPDF(const std::string &filename)
  {    
    FileName = filename;

    // Make a page of 19*15 cm
    PageWidth = cm2pt(19);
    PageHeight = cm2pt(15);
  
    // Axis origo is at 0.1 of PageWidth and PageHeight
    x_origo = PageWidth*0.1;
    y_origo = PageHeight*0.1;

    // Axis cover 80% of page
    AxisWidth = PageWidth*0.8;
    AxisHeight = PageHeight*0.8;

    // Tick Lengths in percentage of axis
    xTickLength = 0.01;
    yTickLength = 0.01;

    AxisLineWidth = 0.5;

    // grid is of by default
    GridLineWidth = 0.1;
    GridLineStyle = DOTTED;
    grid = false;

    // Default font name and size
    dFontName = "Helvetica";
    dFontSize = 10;

    dLineWidth = 0.3;
    dLineStyle = SOLID; // solid line

    MarkerSize = 0.01; // size of makers

    Nplots = 0;

    // Clear all labels and title
    XLabel = "";
    YLabel = "";
    Title = "";

  }



  std::string to_str_fixed(const double &input_value)
  {
    int prec = -1;
    double value;
    std::string output;

    do
      {
	prec++;

	// convert double to string with precision prec
	std::ostringstream out;
	out.precision(prec);
	out.setf(std::ostringstream::fixed, std::ostringstream::floatfield);
	out << input_value;
	output = out.str();

	// convert back to double
	std::istringstream in(output);
	in.setf(std::istringstream::fixed, std::istringstream::floatfield);
	in >> value;
      
      } while ( std::abs(value - input_value) > 5*eps );

    return output;
  }
  


  void PlotPDF::set_axis(const double xmin, const double xmax, const double ymin, const double ymax)
  {
    double range, expbase, tickspace;
    int i;
    vec a = "1 2 5 10";

    XMin = xmin; XMax = xmax; YMin = ymin; YMax = ymax;

    dx = AxisWidth/(XMax-XMin);
    dy = AxisHeight/(YMax-YMin);

    // -------------- x-axis ---------------------------------------------
    // Find x-tick step size
    range = (XMax-XMin);
    expbase = itpp::pow10( floor( log10(range/10) ) );
    i = 0; while (a(i)*expbase < range/10) { i++; }
    tickspace = a(i)*expbase;
  
    int Nx = floor_i( (XMax - round_i(XMin/tickspace)*tickspace)/tickspace + eps ) + 1;

    x_ticks.set_size(Nx, false);
    for (i=0; i<Nx; i++)
      x_ticks(i) = round_i(XMin/tickspace)*tickspace + i*tickspace;
  
    // labels
    xtick_labels.set_size(x_ticks.size(), false);
    for (i=0; i<x_ticks.size(); i++) {
      xtick_labels(i) = to_str_fixed(x_ticks(i));
    }

    //------------- y-axis -------------------------------------------------
    // Find y-tick step size
    range = (YMax-YMin);
    expbase = itpp::pow10( floor( log10(range/10) ) );
    i = 0; while (a(i)*expbase < range/10) { i++; }
    tickspace = a(i)*expbase;
  
    int Ny = floor_i( (YMax - round_i(YMin/tickspace)*tickspace)/tickspace + eps) + 1;

    y_ticks.set_size(Ny, false);
    for (i=0; i<Ny; i++)
      y_ticks(i) = round_i(YMin/tickspace)*tickspace + i*tickspace;
 
    // labels
    ytick_labels.set_size(y_ticks.size(), false);
    for (int i=0; i<y_ticks.size(); i++) {
      ytick_labels(i) = to_str_fixed(y_ticks(i));
    }
  }

  bool PlotPDF::in_range(const double x, const double y)
  {
    return ( ( (x >= XMin) && (x <= XMax) ) && ( (y >= YMin) && (y <= YMax) ) );
  }


  void PlotPDF::draw_line(const double x0, const double y0, const double x1, const double y1)
  {
    canvas->MoveTo( x2pt(x0), y2pt(y0) );
    canvas->LineTo( x2pt(x1), y2pt(y1) );
    canvas->Stroke();
  }

  void PlotPDF::move_to(const double x, const double y)
  {
    canvas->MoveTo( x2pt(x), y2pt(y) );
  }

  void PlotPDF::line_to(const double x, const double y)
  {
    canvas->LineTo( x2pt(x), y2pt(y) );
  }

  void PlotPDF::stroke_line()
  {
    canvas->Stroke();
  }


  void PlotPDF::draw_marker(const double x, const double y, const MarkerType marker)
  {
    double xwidth = MarkerSize*(XMax-XMin);
    double ywidth = MarkerSize*(YMax-YMin);

    switch (marker) {
    case SQUARE: // square marker
      canvas->Rectangle(x2pt(x-xwidth/2), y2pt(y-ywidth/2), xwidth*dx, ywidth*dy);
      canvas->Stroke();
      break;

    default:
      break;
    };

  }

  void PlotPDF::set_font(const std::string &fontname, int fontsize)
  {
    canvas->SetFontAndSize(fontname.c_str(), fontsize);
  }


  void PlotPDF::draw_text(const double x, const double y, const std::string &text)
  {
    canvas->BeginText();
    canvas->MoveTextPos(x2pt(x), y2pt(y));
    canvas->ShowText(text.c_str());
    canvas->EndText();
  }

  void PlotPDF::draw_text_centered(const double x, const double y, const std::string &text)
  {
    double width = canvas->TextWidth(text.c_str());
    canvas->BeginText();
    canvas->MoveTextPos(x2pt(x)-width/2, y2pt(y));
    canvas->ShowText(text.c_str());
    canvas->EndText();
  }

  void PlotPDF::draw_text_right_justified(const double x, const double y, const std::string &text)
  {
    double width = canvas->TextWidth(text.c_str());
    canvas->BeginText();
    canvas->MoveTextPos(x2pt(x)-width, y2pt(y));
    canvas->ShowText(text.c_str());
    canvas->EndText();
  }

  void PlotPDF::draw_vertical_text_centered(const double x, const double y, const std::string &text)
  {
    double width = canvas->TextWidth(text.c_str());
    canvas->BeginText();
    canvas->SetTextMatrix(0.0, 1.0, -1.0, 0.0, x2pt(x), y2pt(y)-width/2);
    canvas->ShowText(text.c_str());
    canvas->EndText();
  }

  void PlotPDF::draw_title()
  {
    if (Title.size() > 0)
      draw_text_centered( (XMax+XMin)/2, YMax + dFontSize*1.0/dy, Title);
  }

  void PlotPDF::draw_xlabel()
  {
    if (XLabel.size() > 0)
      draw_text_centered( (XMax+XMin)/2, YMin - dFontSize*2.6/dy, XLabel);
  }

  void PlotPDF::draw_ylabel()
  {
    if (YLabel.size() > 0)
      draw_vertical_text_centered( XMin - dFontSize*2.6/dx, (YMin+YMax)/2, YLabel);
  }
 

  // Draw x-ticks
  void PlotPDF::draw_xticks()
  {
    double xtick_length = xTickLength*(XMax-XMin);
    set_font(dFontName.c_str(), dFontSize);


    // position in y coordinates of text
    double y_pos = YMin - dFontSize*1.2/dy;
    // How much to move text to left

    for (int i=0; i<x_ticks.size(); i++) {
      draw_line(x_ticks(i), YMin, x_ticks(i), YMin+xtick_length);
      draw_line(x_ticks(i), YMax-xtick_length, x_ticks(i), YMax);
      draw_text_centered(x_ticks(i), y_pos, xtick_labels(i));
    }
  }

  // Draw y-ticks
  void PlotPDF::draw_yticks()
  {
    double ytick_length = yTickLength*(YMax-YMin);
    set_font(dFontName.c_str(), dFontSize);

    // position in x coordinates of text
    double x_pos = XMin - dFontSize*0.4/dx;
    // How much to move text down
    double y_delta = dFontSize*0.4/dy;

    for (int i=0; i<y_ticks.size(); i++) {
      draw_line(XMin, y_ticks(i), XMin+ytick_length, y_ticks(i));
      draw_line(XMax-ytick_length, y_ticks(i), XMax, y_ticks(i));
      draw_text_right_justified(x_pos, y_ticks(i)-y_delta, ytick_labels(i));
    }
  }

  // Draw grid
  void PlotPDF::draw_grid()
  {
    if (grid == true) {
      double xtick_length = xTickLength*(XMax-XMin);
      double ytick_length = yTickLength*(YMax-YMin);

      draw_line_width(GridLineWidth);
      draw_line_style(GridLineStyle);

      for (int i=0; i<x_ticks.size(); i++) {
	draw_line(x_ticks(i), YMin, x_ticks(i), YMax-xtick_length);
      }

      for (int i=0; i<y_ticks.size(); i++) {
	draw_line(XMin, y_ticks(i), XMax-ytick_length, y_ticks(i));
      }
    }
  }


  int PlotPDF::set_plot(const vec &x, const vec &y)
  {
    Nplots++;
    x_values.set_size(Nplots, true);
    y_values.set_size(Nplots, true);
    line_width.set_size(Nplots, true);
    line_style.set_size(Nplots, true);
    line_marker.set_size(Nplots, true);

    x_values(Nplots-1) = x;
    y_values(Nplots-1) = y;

    line_width(Nplots-1) = dLineWidth;
    line_style(Nplots-1) = dLineStyle;
    line_marker(Nplots-1) = NONE;

    return Nplots;
  }

  void PlotPDF::set_line_style(const int plot_number, const LineStyle style)
  {
    // check for valid style
    line_style(plot_number) = style;
  }

  void PlotPDF::set_line_marker(const int plot_number, const MarkerType marker)
  {
    // check for valid style
    line_marker(plot_number) = marker;
  }

  void PlotPDF::set_line_width(const int plot_number, const double width)
  {
    // check for valid width
    line_width(plot_number) = width;
  }
  
  void PlotPDF::draw_line_style(const LineStyle style)
  {
    // check for valid style
    switch (style) {
    case SOLID: // solid line
      canvas->SetDash(0,0,0);
      break;

    case DASHED: // dashed line
      canvas->SetDash(6,6,0);
      break;
      
    case DOTTED: // dotted line
      canvas->SetDash(1,6,0);
      break;

    default: // solid line
      canvas->SetDash(0,0,0);
      break;
    };
  }


  void PlotPDF::draw_line_width(const double width)
  {
    // check for valid width
    canvas->SetLineWidth(width);
  }
  

  // make a funtion add_point(x,y) that stores points on a line segment temporary
  // then plot all these and possible the markers
  // make this a list of points to step through
  void PlotPDF::draw_all_plots()
  {
    int i, j;

    for (i=0; i<Nplots; i++) {
      
      draw_line_style(line_style(i));
      draw_line_width(line_width(i));
      
      j = 0;
      // evaluate all points in vector
      while (j<x_values(i).size()) {
	// Find first point inside axis
	while ( j<x_values(i).size() && !in_range(x_values(i)(j), y_values(i)(j)) ) { j++; }

      	// if first point is out of range create a virtual point on axis
	if( j<x_values(i).size() && (j>0) && !in_range(x_values(i)(j-1), y_values(i)(j-1)) ) {

	  double x0 = x_values(i)(j-1), x1 = x_values(i)(j);
	  double y0 = y_values(i)(j-1), y1 = y_values(i)(j);

	  // equation of line y = kx + m
	  double k = (y1-y0) / (x1-x0); // slope
	  double m = y1-k*x1; // (x1,y1) is inside

	  if( in_range( (YMax-m)/k, YMax) && (y0>=YMax) )
	    move_to((YMax-m)/k, YMax);
	  else if( in_range( (YMin-m)/k, YMin) && (y0<=YMin))
	    move_to((YMin-m)/k, YMin);
	  else if( in_range( XMin, k*XMin+m) && (x0<=XMin) )
	    move_to(XMax, k*XMax+m);

	} else if (j<x_values(i).size())
	  move_to( x_values(i)(j), y_values(i)(j) );
	else
	  break; // did not find any point inside range
      
	// Draw line until point is outside range
	j++;
	while ( j<x_values(i).size() && in_range(x_values(i)(j), y_values(i)(j)) )
	  { line_to( x_values(i)(j), y_values(i)(j) ); j++; }

	// if last point is out of range create a virtual point on axis
	if( j<x_values(i).size() && (j>0) && !in_range(x_values(i)(j), y_values(i)(j)) ) {

	  double x0 = x_values(i)(j-1), x1 = x_values(i)(j);
	  double y0 = y_values(i)(j-1), y1 = y_values(i)(j);

	  // equation of line y = kx + m
	  double k = (y1-y0) / (x1-x0); // slope
	  double m = y0-k*x0; // (x0,y0) is inside

	  if( in_range( (YMax-m)/k, YMax) && (y1>=YMax) )
	    line_to((YMax-m)/k, YMax);
	  else if( in_range( (YMin-m)/k, YMin) && (y1<=YMin))
	    line_to((YMin-m)/k, YMin);
	  else if( in_range( XMax, k*XMax+m) && (x1>=XMax) )
	    line_to(XMax, k*XMax+m);
	}
      
	stroke_line();


      } // while j

    } // for i
  }


  void PlotPDF::save_pdf(void)
  {
    std::string filename = FileName+".pdf"; // add correct suffix

    /* Start creating PDF document. */
    doc = new PdfDoc();

    try {
      /* Add a new page object. */
      doc->NewDoc();

      /* Add Helvetica Font. */
      doc->AddType1Font(new PdfHelveticaFontDef());
      doc->AddType1Font(new PdfHelveticaObliqueFontDef());
    
      page = doc->AddPage();
      page->SetMediaBox(PdfBox(0, 0, (int)PageWidth, (int)PageHeight));
      canvas = page->Canvas();

    
      /* Axis box. */
      canvas->SetLineWidth(AxisLineWidth);
      canvas->Rectangle(x_origo, y_origo, AxisWidth, AxisHeight);
      canvas->Stroke();

      // Draw axis ticks
      draw_xticks();
      draw_yticks();
      
      // Draw grid
      draw_grid();

      draw_title();
      draw_xlabel();
      draw_ylabel();

      // Draw all plots
      canvas->SetLineWidth(dLineWidth);
      draw_all_plots();


      /* Save the document to a file */
      doc->WriteToFile(filename.c_str());
      doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
      fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
  }


  // save in generic format
  void PlotPDF::save()
  {
    std::string filename = FileName + ".itp";
    it_file file(filename);

    file << Name("PageWidth") << PageWidth;
    file << Name("PageHeight") << PageHeight;
    file << Name("AxisWidth") << AxisWidth;
    file << Name("AxisHeight") << AxisHeight;
    file << Name("xTickLength") << xTickLength;
    file << Name("yTickLength") << yTickLength;
    file << Name("AxisLineWidth") << AxisLineWidth;
    file << Name("x_origo") << x_origo;
    file << Name("y_origo") << y_origo;
    file << Name("GridLineWidth") << GridLineWidth;
    file << Name("GridLineStyle") << GridLineStyle;


    file << Name("dFontName") << dFontName;
    file << Name("dFontSize") << dFontSize;
    file << Name("dLineWidth") << dLineWidth;

    // dLineStyle ???
    file << Name("MarkerSize") << MarkerSize;
    
    file << Name("x_ticks") << x_ticks;
    file << Name("y_ticks") << y_ticks;
    file << Name("xtick_labels") << xtick_labels;
    file << Name("ytick_labels") << ytick_labels;
    
    file << Name("XMin") << XMin;
    file << Name("XMax") << XMax;
    file << Name("YMin") << YMin;
    file << Name("YMax") << YMax;
    file << Name("dx") << dx;
    file << Name("dy") << dy;
    file << Name("grid") << grid;

    file << Name("XLabel") << XLabel;
    file << Name("YLabel") << YLabel;
    file << Name("Title") << Title;

    file << Name("Nplots") << Nplots;
    file << Name("x_values") << x_values;
    file << Name("y_values") << y_values;
    file << Name("line_width") << line_width;
    
    // line_style ???
    // line_marker ???

    file << Name("line_color") << line_color;
  }

  // save in generic format
  void PlotPDF::load()
  {
    std::string filename = FileName + ".itp";
    it_file file(filename);

    file >> Name("PageWidth") >> PageWidth;
    file >> Name("PageHeight") >> PageHeight;
    file >> Name("AxisWidth") >> AxisWidth;
    file >> Name("AxisHeight") >> AxisHeight;
    file >> Name("xTickLength") >> xTickLength;
    file >> Name("yTickLength") >> yTickLength;
    file >> Name("AxisLineWidth") >> AxisLineWidth;
    file >> Name("x_origo") >> x_origo;
    file >> Name("y_origo") >> y_origo;
    file >> Name("GridLineWidth") >> GridLineWidth;
    
    // GridLineStyle ???


    file >> Name("dFontName") >> dFontName;
    file >> Name("dFontSize") >> dFontSize;
    file >> Name("dLineWidth") >> dLineWidth;

    // dLineStyle ???
    file >> Name("MarkerSize") >> MarkerSize;
    
    file >> Name("x_ticks") >> x_ticks;
    file >> Name("y_ticks") >> y_ticks;
    file >> Name("xtick_labels") >> xtick_labels;
    file >> Name("ytick_labels") >> ytick_labels;
    
    file >> Name("XMin") >> XMin;
    file >> Name("XMax") >> XMax;
    file >> Name("YMin") >> YMin;
    file >> Name("YMax") >> YMax;
    file >> Name("dx") >> dx;
    file >> Name("dy") >> dy;

    // grid ???

    file >> Name("XLabel") >> XLabel;
    file >> Name("YLabel") >> YLabel;
    file >> Name("Title") >> Title;

    file >> Name("Nplots") >> Nplots;
    file >> Name("x_values") >> x_values;
    file >> Name("y_values") >> y_values;
    file >> Name("line_width") >> line_width;
    
    // line_style ???
    // line_marker ???

    file >> Name("line_color") >> line_color;
  }

} //namespace itpp

#endif // HAVE_HARU
