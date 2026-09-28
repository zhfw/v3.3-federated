import streamlit as st
import pandas as pd
from pathlib import Path
from predictor_core import predict
from report_generator import build_pdf

st.set_page_config(page_title='CFPS Manufacturability Predictor v3.3 Integrated',page_icon='🧬',layout='wide')
st.markdown('''<style>
.block-container{padding-top:1.2rem;max-width:1380px}.badge{display:inline-block;padding:5px 10px;border-radius:999px;font-weight:700}.go{background:#dcfae6;color:#067647}.hold{background:#fee4e2;color:#b42318}.prov{background:#fef0c7;color:#b54708}.warnbox{border-left:5px solid #d92d20;background:#fff5f4;padding:12px;border-radius:8px}.muted{color:#667085;font-size:.88rem}.card{border:1px solid #e4e7ec;border-radius:12px;padding:12px;background:#fff}
</style>''',unsafe_allow_html=True)
ROOT=Path(__file__).resolve().parent
DEMO_PROTEIN=(ROOT/'data'/'demo_protein.txt').read_text().strip();DEMO_DNA=(ROOT/'data'/'demo_dna.txt').read_text().strip()

st.title('CFPS Manufacturability Predictor v3.3')
st.caption('v3.3 trained ML backend + v2.5 warning/release/PDF layer + federated-ready architecture')

tabs=st.tabs(['Predictor','Warnings & solutions','Construct diagnostics','PDF report','Federated network'])
with tabs[0]:
    left,right=st.columns([1.12,.88])
    with left:
        protein=st.text_area('Protein amino-acid sequence',DEMO_PROTEIN,height=220)
        dna=st.text_area('Coding DNA (optional)',DEMO_DNA,height=120)
        c1,c2,c3=st.columns(3);system=c1.selectbox('Expression context',['E. coli CFPS / PURE-like','E. coli recombinant production']);temp=c2.selectbox('Planned temperature',['Not specified','16 C','20 C','25 C','30 C','37 C']);chap=c3.selectbox('Chaperone condition',['Not specified','None','Trigger Factor','GroEL/ES','DnaK/J/GrpE'])
        construct_id=st.text_input('Construct ID (optional)','Demo-01')
        if st.button('Run integrated v3.3 prediction',type='primary',use_container_width=True):
            try:
                st.session_state['result']=predict(protein,dna or None);st.session_state['inputs']={'protein':protein,'dna':dna,'system':system,'temperature':temp,'chaperone':chap,'construct_id':construct_id}
            except Exception as e: st.error(str(e))
    r=st.session_state.get('result')
    with right:
        if r:
            rel=r['release'];cls='go' if rel['status']=='GO' else ('hold' if rel['status']=='HOLD' else 'prov');st.markdown(f"<span class='badge {cls}'>{rel['status']}</span>",unsafe_allow_html=True)
            a,b=st.columns(2);a.metric('v2.5 compatibility score','N/A' if rel['combined_score'] is None else f"{rel['combined_score']:.1f}/100");b.metric('v3.3 prioritization',f"{100*r['overall_score']:.1f}/100")
            a,b=st.columns(2);a.metric('Protein gate',f"{rel['protein_score']:.1f}/100",'PASS' if rel['protein_pass'] else 'FAIL');b.metric('DNA gate','N/A' if rel['dna_score'] is None else f"{rel['dna_score']:.1f}/100",'DNA required' if rel['dna_pass'] is None else ('PASS' if rel['dna_pass'] else 'FAIL'))
            if rel['red_warning']: st.markdown("<div class='warnbox'><b>Release warning:</b> combined score is below 85/100. Optimization or additional validation is recommended before release.</div>",unsafe_allow_html=True)
            st.divider();a,b=st.columns(2);a.metric('CFPS soluble fraction',f"{100*r['cfps_soluble_fraction']:.1f}%");b.metric('Aggregation / insoluble',f"{100*r['cfps_aggregation_fraction']:.1f}%")
            a,b=st.columns(2);a.metric('RP3 production',f"{100*r['rp3_production_probability_protein']:.1f}%");b.metric('General solubility',f"{100*r['general_solubility_probability']:.1f}%")
            if r['rp3_production_probability_with_dna'] is not None: st.metric('RP3 production - protein + DNA',f"{100*r['rp3_production_probability_with_dna']:.1f}%")
            a,b=st.columns(2);a.metric('NESG expression',f"{r['nesg_expression_score_0_5']:.2f}/5");b.metric('NESG solubility',f"{r['nesg_solubility_score_0_5']:.2f}/5")
            st.caption(f"Confidence: {r['confidence']} | OOD index: {r['ood_index']:.2f}")
        else: st.info('Run the prediction. A matched demonstration protein/DNA pair is preloaded.')

with tabs[1]:
    r=st.session_state.get('result')
    if not r: st.info('Run a prediction first.')
    else:
        st.subheader('Warnings, consequences, and possible solutions')
        severity_order={'high':0,'medium':1,'low':2}
        for d in sorted(r['diagnostics'],key=lambda x:severity_order.get(x['level'],9)):
            icon={'high':'🔴','medium':'🟠','low':'🟢'}[d['level']]
            with st.expander(f"{icon} {d['title']} - {d['level'].upper()}",expanded=d['level']=='high'):
                st.markdown(f"**Evidence:** {d['detail']}")
                st.markdown(f"**Likely consequence:** {d['consequence']}")
                st.markdown(f"**Possible solution:** {d['solution']}")
                st.markdown(f"**Suggested validation:** {d['validation']}")
        st.subheader('Prioritized experimental validation plan')
        for i,x in enumerate(r['validation_plan'],1): st.write(f'{i}. {x}')
        st.info('Solutions are decision-support suggestions. Sequence redesign should preserve functional/structural constraints and be experimentally validated.')

with tabs[2]:
    r=st.session_state.get('result')
    if not r: st.info('Run a prediction first.')
    else:
        f=r['protein_features'];c1,c2,c3,c4=st.columns(4);c1.metric('Length',int(f['length']));c2.metric('Approx. MW',f"{f['mw_approx']/1000:.1f} kDa");c3.metric('GRAVY',f"{f['gravy']:.3f}");c4.metric('Cysteines',int(f['cys_count']))
        c1,c2,c3,c4=st.columns(4);c1.metric('Hydrophobic fraction',f"{f['hydrophobic_frac']:.1%}");c2.metric('Longest hydro. run',int(f['longest_hydrophobic_run']));c3.metric('Charge proxy',f"{f['charge_proxy']:+.3f}");c4.metric('Entropy',f"{f['entropy']:.2f}")
        with st.expander('All engineered protein features'): st.dataframe(pd.DataFrame({'feature':list(f.keys()),'value':list(f.values())}),hide_index=True,use_container_width=True)
        if r['dna_features']:
            d=r['dna_features'];st.subheader('DNA diagnostics');c1,c2,c3,c4=st.columns(4);c1.metric('GC',f"{d['gc']:.1%}");c2.metric('GC3',f"{d['gc3']:.1%}");c3.metric('Longest homopolymer',int(d['longest_homopolymer']));c4.metric('Start ATG','Yes' if d['start_atg'] else 'No')

with tabs[3]:
    r=st.session_state.get('result'); inp=st.session_state.get('inputs')
    if not r: st.info('Run a prediction first.')
    else:
        st.subheader('PDF report')
        st.write('The report contains the v3.3 predictions, v2.5 release gates, red warnings, likely consequences, possible solutions, validation suggestions, construct diagnostics, and model limitations.')
        pdf=build_pdf(r,inp['protein'],inp['dna'] or None,{'Construct ID':inp.get('construct_id',''),'Expression context':inp.get('system',''),'Planned temperature':inp.get('temperature',''),'Chaperone condition':inp.get('chaperone','')})
        filename=(inp.get('construct_id') or 'construct').replace(' ','_')+'_v33_report.pdf'
        st.download_button('Download PDF report',pdf,file_name=filename,mime='application/pdf',type='primary',use_container_width=True)
        st.caption('Release rule shown in PDF: 70% protein + 30% DNA; GO only if combined >=85, protein >=72, DNA >=78.')

with tabs[4]:
    st.subheader('Three-party CPU federated mode')
    st.code('Party A: private sequence -> frozen RP3Net -> 256-D STP -> local head  \\nParty B: private sequence -> frozen RP3Net -> 256-D STP -> local head   > FedAvg / FedProx -> global v3.3 head\nParty C: private sequence -> frozen RP3Net -> 256-D STP -> local head  /')
    m1,m2,m3,m4=st.columns(4);m1.metric('RP3 FedAvg AUROC','0.773');m2.metric('RP3 AUPRC','0.396');m3.metric('CFPS FedProx R2','0.402');m4.metric('CFPS MAE','0.190')
    st.markdown('**Privacy boundary:** raw sequences, coding DNA, experimental measurements, and local STP feature stores remain local. Only agreed model-head updates and aggregate round metadata are exchanged.')
    st.info('This page is a local simulator/prototype. Production federation should use authenticated transport and secure aggregation rather than uploading confidential site datasets through Streamlit.')
