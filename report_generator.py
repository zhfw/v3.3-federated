from io import BytesIO
from datetime import datetime
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib.units import mm

SEVERITY_COLORS = {
    'high': colors.HexColor('#FDE8E7'),
    'medium': colors.HexColor('#FFF4D6'),
    'low': colors.HexColor('#E7F6EC'),
}

def _pct(x):
    return 'N/A' if x is None else f'{100*x:.1f}%'

def _score(x):
    return 'N/A' if x is None else f'{x:.1f}/100'

def build_pdf(result, protein_seq, dna_seq=None, metadata=None):
    metadata = metadata or {}
    bio = BytesIO()
    doc = SimpleDocTemplate(
        bio, pagesize=A4, rightMargin=15*mm, leftMargin=15*mm,
        topMargin=15*mm, bottomMargin=15*mm,
        title='CFPS Manufacturability Predictor v3.3 + v2.5 Report'
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='TitleCenter', parent=styles['Title'], alignment=TA_CENTER, fontSize=19, leading=23, spaceAfter=8))
    styles.add(ParagraphStyle(name='Subtle', parent=styles['BodyText'], textColor=colors.HexColor('#667085'), fontSize=8.5, leading=11))
    styles.add(ParagraphStyle(name='H2x', parent=styles['Heading2'], fontSize=13, leading=16, spaceBefore=8, spaceAfter=5))
    styles.add(ParagraphStyle(name='Small', parent=styles['BodyText'], fontSize=8.5, leading=11))
    styles.add(ParagraphStyle(name='Warn', parent=styles['BodyText'], fontSize=9, leading=12, textColor=colors.HexColor('#B42318')))

    story=[]
    story.append(Paragraph('CFPS Manufacturability Predictor', styles['TitleCenter']))
    story.append(Paragraph('v3.3 trained prediction backend + v2.5 diagnostic/release layer', styles['Subtle']))
    story.append(Spacer(1, 6))

    release=result['release']
    status=release['status']
    status_color={'GO':'#067647','HOLD':'#B42318','PROVISIONAL':'#B54708'}.get(status,'#344054')
    summary_data=[
        ['Release status', f"<b><font color='{status_color}'>{status}</font></b>"],
        ['v2.5 compatibility score', _score(release['combined_score'])],
        ['Protein gate score', _score(release['protein_score'])],
        ['DNA gate score', _score(release['dna_score'])],
        ['v3.3 overall prioritization', f"{100*result['overall_score']:.1f}/100"],
        ['Confidence', result['confidence']],
    ]
    tbl=Table([[Paragraph(str(a),styles['Small']),Paragraph(str(b),styles['Small'])] for a,b in summary_data], colWidths=[62*mm,105*mm])
    tbl.setStyle(TableStyle([
        ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#F2F4F7')),
        ('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#D0D5DD')),
        ('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)
    ]))
    story.append(tbl)
    story.append(Spacer(1,8))
    if release['red_warning']:
        story.append(Paragraph('<b>WARNING:</b> The combined release score is below 85/100. Treat this construct as requiring optimization or additional validation before release.', styles['Warn']))
        story.append(Spacer(1,6))

    story.append(Paragraph('1. Prediction summary', styles['H2x']))
    preds=[
        ['Endpoint','Prediction'],
        ['CFPS soluble fraction', _pct(result['cfps_soluble_fraction'])],
        ['Operational aggregation / insoluble fraction', _pct(result['cfps_aggregation_fraction'])],
        ['General solubility probability', _pct(result['general_solubility_probability'])],
        ['RP3 production probability - protein', _pct(result['rp3_production_probability_protein'])],
        ['RP3 production probability - protein + DNA', _pct(result['rp3_production_probability_with_dna'])],
        ['NESG expression score', f"{result['nesg_expression_score_0_5']:.2f}/5"],
        ['NESG solubility score', f"{result['nesg_solubility_score_0_5']:.2f}/5"],
        ['OOD index', f"{result['ood_index']:.2f}"],
    ]
    t=Table([[Paragraph(str(x),styles['Small']) for x in row] for row in preds], colWidths=[105*mm,62*mm])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAECF0')),('GRID',(0,0),(-1,-1),0.35,colors.HexColor('#D0D5DD')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
    story.append(t)

    story.append(Paragraph('2. v2.5 compatibility/release logic', styles['H2x']))
    story.append(Paragraph('The compatibility layer preserves the earlier v2.5 release rule: combined score = 70% protein score + 30% DNA score. GO requires combined >= 85, protein >= 72, and DNA >= 78. When DNA is absent, the result is PROVISIONAL and cannot receive final GO.', styles['Small']))
    gates=[['Gate','Value','Requirement','Pass'],
           ['Combined',_score(release['combined_score']),'>=85',str(release['combined_pass'])],
           ['Protein',_score(release['protein_score']),'>=72',str(release['protein_pass'])],
           ['DNA',_score(release['dna_score']),'>=78',str(release['dna_pass']) if release['dna_score'] is not None else 'DNA missing']]
    t=Table([[Paragraph(str(x),styles['Small']) for x in row] for row in gates],colWidths=[38*mm,42*mm,42*mm,45*mm])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAECF0')),('GRID',(0,0),(-1,-1),0.35,colors.HexColor('#D0D5DD')),('VALIGN',(0,0),(-1,-1),'TOP')]))
    story.append(Spacer(1,4)); story.append(t)

    story.append(Paragraph('3. Warnings, likely consequences, and corrective actions', styles['H2x']))
    warning_rows=[[Paragraph('<b>Severity</b>',styles['Small']),Paragraph('<b>Issue</b>',styles['Small']),Paragraph('<b>Evidence / consequence</b>',styles['Small']),Paragraph('<b>Recommended action</b>',styles['Small'])]]
    for w in result['diagnostics']:
        warning_rows.append([
            Paragraph(w['level'].upper(),styles['Small']), Paragraph(w['title'],styles['Small']),
            Paragraph(w['detail']+' '+w.get('consequence',''),styles['Small']), Paragraph(w['solution'],styles['Small'])
        ])
    wt=Table(warning_rows,colWidths=[22*mm,42*mm,62*mm,50*mm],repeatRows=1)
    style=[('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAECF0')),('GRID',(0,0),(-1,-1),0.3,colors.HexColor('#D0D5DD')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]
    for i,w in enumerate(result['diagnostics'], start=1): style.append(('BACKGROUND',(0,i),(0,i),SEVERITY_COLORS.get(w['level'],colors.white)))
    wt.setStyle(TableStyle(style)); story.append(wt)

    story.append(Paragraph('4. Experimental follow-up suggestions', styles['H2x']))
    for x in result['validation_plan']:
        story.append(Paragraph('- '+x, styles['Small']))

    story.append(Paragraph('5. Construct diagnostics', styles['H2x']))
    f=result['protein_features']
    diag=[['Feature','Value'],['Length',str(int(f['length']))],['Approx. molecular weight',f"{f['mw_approx']/1000:.1f} kDa"],['GRAVY',f"{f['gravy']:.3f}"],['Hydrophobic fraction',f"{f['hydrophobic_frac']:.1%}"],['Longest hydrophobic run',str(int(f['longest_hydrophobic_run']))],['Charge proxy',f"{f['charge_proxy']:+.3f}"],['Cysteines',str(int(f['cys_count']))],['Sequence entropy',f"{f['entropy']:.2f}"]]
    if result['dna_features']:
        d=result['dna_features']; diag += [['DNA GC',f"{d['gc']:.1%}"],['DNA GC3',f"{d['gc3']:.1%}"],['Longest DNA homopolymer',str(int(d['longest_homopolymer']))]]
    dt=Table([[Paragraph(str(x),styles['Small']) for x in row] for row in diag],colWidths=[90*mm,77*mm])
    dt.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAECF0')),('GRID',(0,0),(-1,-1),0.3,colors.HexColor('#D0D5DD'))]));story.append(dt)

    story.append(PageBreak())
    story.append(Paragraph('6. Input record and model notes', styles['H2x']))
    story.append(Paragraph(f"Protein sequence length: {len(''.join(protein_seq.split()))} aa",styles['Small']))
    story.append(Paragraph(f"Coding DNA supplied: {'yes' if dna_seq else 'no'}",styles['Small']))
    if metadata:
        for k,v in metadata.items(): story.append(Paragraph(f"{k}: {v}",styles['Small']))
    story.append(Spacer(1,6))
    story.append(Paragraph('Model limitations', styles['H2x']))
    for x in result['limitations']:
        story.append(Paragraph('- '+x,styles['Small']))
    story.append(Spacer(1,8))
    story.append(Paragraph('Interpretation note: model outputs are prioritization and risk-estimation aids. They do not replace experimental validation. Absolute CFPS yield and CHO secretion/titer are not claimed unless a trained quantitative head is available.', styles['Subtle']))
    story.append(Spacer(1,8))
    story.append(Paragraph('Generated: '+datetime.now().strftime('%Y-%m-%d %H:%M'),styles['Subtle']))

    doc.build(story)
    return bio.getvalue()
