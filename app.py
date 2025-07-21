# Simplified Protein Structure-Function Web Application
# Compatible with Python 3.13

from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
import os
import json
import numpy as np
from datetime import datetime
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# Configuration
UPLOAD_FOLDER = 'uploads'
MODELS_FOLDER = 'models'
RESULTS_FOLDER = 'results'

# Ensure directories exist
for folder in [UPLOAD_FOLDER, MODELS_FOLDER, RESULTS_FOLDER]:
    os.makedirs(folder, exist_ok=True)

class SimpleProteinAnalyzer:
    """
    Simplified protein analyzer that works without complex dependencies
    """
    
    def __init__(self):
        self.dssp_states = ['H', 'G', 'I', 'E', 'B', 'T', 'S', 'C']  # 8-state DSSP
        self.state_descriptions = {
            'H': 'α-helix (Alpha helix)',
            'G': '3₁₀-helix (3-10 helix)', 
            'I': 'π-helix (Pi helix)',
            'E': 'β-strand (Extended strand)',
            'B': 'β-bridge (Beta bridge)',
            'T': 'Turn',
            'S': 'Bend', 
            'C': 'Coil/Loop'
        }
        self.functional_implications = {
            'H': 'Often involved in protein-protein interactions, binding sites, and structural stability',
            'G': 'Tight turns, often at protein surfaces, involved in ligand binding',
            'I': 'Rare, usually in membrane proteins or special structural contexts',
            'E': 'Critical for β-sheets, enzyme active sites, protein-protein interfaces',
            'B': 'Bridge between β-strands, important for sheet topology',
            'T': 'Flexible regions, often involved in conformational changes',
            'S': 'Flexible bends, important for protein dynamics',
            'C': 'Highly flexible, often contains binding sites and regulatory regions'
        }
        
        # Amino acid propensities (Chou-Fasman parameters)
        self.helix_propensity = {
            'A': 1.42, 'E': 1.51, 'L': 1.21, 'M': 1.45, 'Q': 1.11, 'R': 0.98, 'K': 1.16,
            'D': 1.01, 'N': 0.67, 'S': 0.77, 'T': 0.83, 'C': 0.70, 'V': 1.06, 'I': 1.08,
            'F': 1.13, 'Y': 0.69, 'W': 1.08, 'H': 1.00, 'P': 0.57, 'G': 0.57
        }
        
        self.sheet_propensity = {
            'V': 1.70, 'I': 1.60, 'Y': 1.47, 'F': 1.38, 'W': 1.37, 'L': 1.30, 'T': 1.19,
            'C': 1.19, 'A': 0.83, 'R': 0.93, 'G': 0.75, 'M': 1.05, 'N': 0.89, 'P': 0.55,
            'S': 0.75, 'H': 0.87, 'Q': 1.10, 'E': 0.37, 'K': 0.74, 'D': 0.54
        }
        
        self.turn_propensity = {
            'N': 1.56, 'G': 1.56, 'P': 1.52, 'S': 1.43, 'D': 1.46, 'T': 0.96, 'H': 0.95,
            'C': 1.19, 'Y': 1.14, 'K': 1.01, 'Q': 0.98, 'R': 0.95, 'W': 0.96, 'A': 0.66,
            'M': 0.60, 'I': 0.47, 'L': 0.59, 'V': 0.50, 'F': 0.60, 'E': 0.74
        }
    
    def predict_structure_function(self, sequence):
        """
        Predict secondary structure using Chou-Fasman rules and analyze function
        """
        try:
            # Predict secondary structure
            predicted_structure = self._predict_secondary_structure(sequence)
            
            # Calculate confidence based on propensity scores
            confidence = self._calculate_confidence(sequence, predicted_structure)
            
            # Analyze structure-function relationships
            functional_analysis = self._analyze_structure_function(sequence, predicted_structure)
            
            # Create detailed results
            results = {
                'sequence': sequence,
                'structure': predicted_structure,
                'confidence': confidence,
                'functional_analysis': functional_analysis,
                'detailed_prediction': self._create_detailed_prediction(sequence, predicted_structure),
                'statistics': self._calculate_structure_statistics(predicted_structure),
                'visualization_data': self._prepare_visualization_data(sequence, predicted_structure)
            }
            
            return results
            
        except Exception as e:
            logger.error(f"Error in structure-function prediction: {e}")
            return self._fallback_prediction(sequence)
    
    def _predict_secondary_structure(self, sequence):
        """
        Predict secondary structure using Chou-Fasman method
        """
        structure = []
        
        for i, aa in enumerate(sequence):
            # Get propensities for current amino acid
            h_prop = self.helix_propensity.get(aa, 1.0)
            e_prop = self.sheet_propensity.get(aa, 1.0) 
            t_prop = self.turn_propensity.get(aa, 1.0)
            
            # Consider local context (simplified sliding window)
            context_h = self._get_context_propensity(sequence, i, self.helix_propensity)
            context_e = self._get_context_propensity(sequence, i, self.sheet_propensity)
            context_t = self._get_context_propensity(sequence, i, self.turn_propensity)
            
            # Weighted decision
            h_score = h_prop * 0.7 + context_h * 0.3
            e_score = e_prop * 0.7 + context_e * 0.3
            t_score = t_prop * 0.7 + context_t * 0.3
            
            # Apply thresholds and rules
            if h_score > 1.03 and h_score > e_score:
                structure.append('H')
            elif e_score > 1.03 and e_score > h_score:
                structure.append('E')
            elif t_score > 1.0 or aa in ['P', 'G']:
                structure.append('T')
            else:
                structure.append('C')
        
        # Post-process to add other DSSP states
        processed_structure = self._add_dssp_states(structure, sequence)
        
        return ''.join(processed_structure)
    
    def _get_context_propensity(self, sequence, pos, propensity_dict):
        """
        Calculate average propensity in local window
        """
        window_size = 5
        half_window = window_size // 2
        
        start = max(0, pos - half_window)
        end = min(len(sequence), pos + half_window + 1)
        
        propensities = []
        for i in range(start, end):
            if i != pos:  # Exclude current position
                aa = sequence[i]
                prop = propensity_dict.get(aa, 1.0)
                propensities.append(prop)
        
        return sum(propensities) / len(propensities) if propensities else 1.0
    
    def _add_dssp_states(self, structure, sequence):
        """
        Add additional DSSP states (G, I, B, S) based on context
        """
        processed = []
        
        for i, ss in enumerate(structure):
            aa = sequence[i]
            
            if ss == 'H':
                # Check for 3-10 helix (G) - Pro often breaks alpha helix
                if i > 0 and i < len(sequence) - 1:
                    if aa == 'P' or sequence[i-1] == 'P' or sequence[i+1] == 'P':
                        processed.append('G')
                    elif aa in ['F', 'W', 'Y'] and self._is_near_membrane_region(sequence, i):
                        processed.append('I')  # Pi helix in membrane proteins
                    else:
                        processed.append('H')
                else:
                    processed.append('H')
                    
            elif ss == 'E':
                # Check for beta bridge (B) - isolated strands
                is_isolated = self._is_isolated_strand(structure, i)
                if is_isolated:
                    processed.append('B')
                else:
                    processed.append('E')
                    
            elif ss == 'T':
                # Distinguish between turn (T) and bend (S)
                if aa in ['G', 'P'] or self._has_sharp_turn(sequence, i):
                    processed.append('T')
                else:
                    processed.append('S')
            else:
                processed.append('C')
        
        return processed
    
    def _is_near_membrane_region(self, sequence, pos):
        """
        Simple check for hydrophobic regions (potential membrane)
        """
        hydrophobic = set(['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'])
        window_size = 15
        start = max(0, pos - window_size // 2)
        end = min(len(sequence), pos + window_size // 2 + 1)
        
        hydrophobic_count = sum(1 for aa in sequence[start:end] if aa in hydrophobic)
        return hydrophobic_count / (end - start) > 0.6
    
    def _is_isolated_strand(self, structure, pos):
        """
        Check if beta strand is isolated (beta bridge)
        """
        window_size = 6
        start = max(0, pos - window_size)
        end = min(len(structure), pos + window_size + 1)
        
        local_structure = structure[start:end]
        strand_count = local_structure.count('E')
        
        return strand_count <= 2  # Isolated if few strands nearby
    
    def _has_sharp_turn(self, sequence, pos):
        """
        Check for amino acids that promote sharp turns
        """
        turn_promoters = set(['P', 'G', 'N', 'D', 'S'])
        if pos == 0 or pos >= len(sequence) - 1:
            return False
        
        return (sequence[pos] in turn_promoters or 
                sequence[pos-1] in turn_promoters or 
                sequence[pos+1] in turn_promoters)
    
    def _calculate_confidence(self, sequence, structure):
        """
        Calculate prediction confidence based on propensity scores
        """
        total_score = 0
        for i, (aa, ss) in enumerate(zip(sequence, structure)):
            if ss in ['H', 'G', 'I']:
                score = self.helix_propensity.get(aa, 1.0)
            elif ss in ['E', 'B']:
                score = self.sheet_propensity.get(aa, 1.0)
            elif ss in ['T', 'S']:
                score = self.turn_propensity.get(aa, 1.0)
            else:
                score = 1.0 - max(
                    self.helix_propensity.get(aa, 1.0) - 1.0,
                    self.sheet_propensity.get(aa, 1.0) - 1.0,
                    0
                )
            
            total_score += min(score, 2.0)  # Cap extreme values
        
        avg_score = total_score / len(sequence)
        confidence = min(max((avg_score - 0.5) * 100, 45), 95)  # Scale to 45-95%
        
        return round(confidence, 2)
    
    def _analyze_structure_function(self, sequence, structure):
        """
        Analyze structure-function relationships
        """
        analysis = {
            'functional_regions': [],
            'structural_motifs': [],
            'binding_potential': [],
            'flexibility_analysis': {},
            'conservation_importance': {}
        }
        
        # Identify functional regions
        current_region = None
        region_start = 0
        
        for i, ss in enumerate(structure):
            if current_region != ss:
                if current_region and i - region_start >= 3:  # Minimum length
                    region = {
                        'type': current_region,
                        'start': region_start + 1,  # 1-indexed
                        'end': i,
                        'length': i - region_start,
                        'sequence': sequence[region_start:i],
                        'description': self.state_descriptions[current_region],
                        'functional_implications': self.functional_implications[current_region]
                    }
                    analysis['functional_regions'].append(region)
                
                current_region = ss
                region_start = i
        
        # Add final region
        if current_region and len(structure) - region_start >= 3:
            region = {
                'type': current_region,
                'start': region_start + 1,
                'end': len(structure),
                'length': len(structure) - region_start,
                'sequence': sequence[region_start:],
                'description': self.state_descriptions[current_region],
                'functional_implications': self.functional_implications[current_region]
            }
            analysis['functional_regions'].append(region)
        
        # Identify structural motifs
        analysis['structural_motifs'] = self._identify_motifs(sequence, structure)
        
        # Analyze binding potential
        analysis['binding_potential'] = self._analyze_binding_potential(sequence, structure)
        
        # Flexibility analysis
        analysis['flexibility_analysis'] = self._analyze_flexibility(structure)
        
        return analysis
    
    def _identify_motifs(self, sequence, structure):
        """
        Identify common structural motifs
        """
        motifs = []
        
        # Helix-turn-helix motif
        import re
        hth_pattern = r'[HGI]{6,}[TS]{2,6}[HGI]{6,}'
        for match in re.finditer(hth_pattern, structure):
            motifs.append({
                'type': 'Helix-Turn-Helix',
                'start': match.start() + 1,
                'end': match.end(),
                'function': 'DNA binding domain, transcriptional regulation',
                'sequence': sequence[match.start():match.end()]
            })
        
        # Beta hairpin
        hairpin_pattern = r'[EB]{3,}[TS]{2,4}[EB]{3,}'
        for match in re.finditer(hairpin_pattern, structure):
            motifs.append({
                'type': 'Beta Hairpin',
                'start': match.start() + 1,
                'end': match.end(),
                'function': 'Structural stability, protein-protein interactions',
                'sequence': sequence[match.start():match.end()]
            })
        
        return motifs
    
    def _analyze_binding_potential(self, sequence, structure):
        """
        Analyze potential binding sites based on structure
        """
        binding_sites = []
        
        # Look for flexible regions (potential binding sites)
        for i, ss in enumerate(structure):
            if ss in ['C', 'T', 'S'] and i < len(sequence) - 2:
                # Check if surrounded by structured regions
                context_structured = False
                if i > 3 and i < len(structure) - 3:
                    before = structure[i-3:i]
                    after = structure[i+1:i+4]
                    if any(s in 'HGIE' for s in before) and any(s in 'HGIE' for s in after):
                        context_structured = True
                
                if context_structured:
                    binding_sites.append({
                        'position': i + 1,
                        'type': 'Flexible loop',
                        'potential': 'High',
                        'reason': 'Flexible region between structured elements',
                        'amino_acid': sequence[i]
                    })
        
        return binding_sites
    
    def _analyze_flexibility(self, structure):
        """
        Analyze protein flexibility based on secondary structure
        """
        flexibility_scores = []
        flexibility_map = {'H': 0.2, 'G': 0.3, 'I': 0.1, 'E': 0.1, 'B': 0.4, 'T': 0.8, 'S': 0.7, 'C': 0.9}
        
        for ss in structure:
            flexibility_scores.append(flexibility_map[ss])
        
        return {
            'average_flexibility': round(sum(flexibility_scores) / len(flexibility_scores), 3),
            'most_flexible_regions': [i+1 for i, score in enumerate(flexibility_scores) if score > 0.7],
            'most_rigid_regions': [i+1 for i, score in enumerate(flexibility_scores) if score < 0.3],
            'flexibility_profile': flexibility_scores
        }
    
    def _create_detailed_prediction(self, sequence, structure):
        """
        Create position-by-position detailed prediction
        """
        detailed = []
        for i, (aa, ss) in enumerate(zip(sequence, structure)):
            detailed.append({
                'position': i + 1,
                'amino_acid': aa,
                'secondary_structure': ss,
                'description': self.state_descriptions[ss],
                'structural_class': self._get_structural_class(ss)
            })
        return detailed
    
    def _get_structural_class(self, ss):
        """
        Get broader structural classification
        """
        if ss in ['H', 'G', 'I']:
            return 'Helical'
        elif ss in ['E', 'B']:
            return 'Extended'
        else:
            return 'Coil'
    
    def _calculate_structure_statistics(self, structure):
        """
        Calculate secondary structure statistics
        """
        counts = {ss: structure.count(ss) for ss in self.dssp_states}
        total = len(structure)
        
        stats = {
            'composition': {ss: {'count': count, 'percentage': round((count/total)*100, 1)} 
                          for ss, count in counts.items()},
            'summary': {
                'total_length': total,
                'helical_content': round(((counts['H'] + counts['G'] + counts['I'])/total)*100, 1),
                'extended_content': round(((counts['E'] + counts['B'])/total)*100, 1),
                'coil_content': round(((counts['C'] + counts['T'] + counts['S'])/total)*100, 1)
            }
        }
        return stats
    
    def _prepare_visualization_data(self, sequence, structure):
        """
        Prepare data for visualization
        """
        return {
            'sequence_positions': list(range(1, len(sequence) + 1)),
            'amino_acids': list(sequence),
            'secondary_structure': list(structure),
            'structure_colors': [self._get_structure_color(ss) for ss in structure],
            'flexibility_profile': [self._get_flexibility_score(ss) for ss in structure]
        }
    
    def _get_structure_color(self, ss):
        """
        Get color for secondary structure visualization
        """
        colors = {
            'H': '#FF0000',  # Red for helix
            'G': '#FF6666',  # Light red for 3-10 helix
            'I': '#990000',  # Dark red for pi helix
            'E': '#0000FF',  # Blue for sheet
            'B': '#6666FF',  # Light blue for bridge
            'T': '#00FF00',  # Green for turn
            'S': '#66FF66',  # Light green for bend
            'C': '#CCCCCC'   # Gray for coil
        }
        return colors.get(ss, '#000000')
    
    def _get_flexibility_score(self, ss):
        """
        Get flexibility score for visualization
        """
        scores = {'H': 0.2, 'G': 0.3, 'I': 0.1, 'E': 0.1, 'B': 0.4, 'T': 0.8, 'S': 0.7, 'C': 0.9}
        return scores.get(ss, 0.5)
    
    def _fallback_prediction(self, sequence):
        """
        Fallback prediction if main model fails
        """
        # Very simple rule-based prediction
        structure = ''
        for aa in sequence:
            if aa in ['A', 'E', 'L', 'M']:
                structure += 'H'
            elif aa in ['V', 'I', 'F', 'Y']:
                structure += 'E'
            elif aa in ['G', 'P']:
                structure += 'T'
            else:
                structure += 'C'
        
        return {
            'sequence': sequence,
            'structure': structure,
            'confidence': 60.0,
            'functional_analysis': {'functional_regions': [], 'note': 'Simplified prediction'},
            'statistics': self._calculate_structure_statistics(structure)
        }

# Initialize the analyzer
analyzer = SimpleProteinAnalyzer()

@app.route('/')
def index():
    """Main page"""
    return render_template('enhanced_index.html')

@app.route('/api/predict_structure_function', methods=['POST'])
def predict_structure_function():
    """Enhanced structure-function prediction endpoint"""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').strip().upper()
        
        if not sequence:
            return jsonify({'error': 'No sequence provided'}), 400
        
        # Validate sequence
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        clean_sequence = sequence.replace('\n', '').replace(' ', '')
        
        if not all(aa in valid_aa for aa in clean_sequence):
            return jsonify({'error': 'Invalid amino acid sequence'}), 400
        
        if len(clean_sequence) < 3:
            return jsonify({'error': 'Sequence too short (minimum 3 amino acids)'}), 400
        
        if len(clean_sequence) > 1000:
            return jsonify({'error': 'Sequence too long (maximum 1000 amino acids)'}), 400
        
        # Perform prediction
        results = analyzer.predict_structure_function(clean_sequence)
        
        # Add metadata
        results['analysis_type'] = 'Chou-Fasman based 8-state prediction'
        results['timestamp'] = datetime.now().isoformat()
        results['sequence_length'] = len(clean_sequence)
        
        return jsonify(results)
        
    except Exception as e:
        logger.error(f"Structure-function prediction error: {e}")
        return jsonify({'error': f'Analysis failed: {str(e)}'}), 500

@app.route('/api/analyze_batch', methods=['POST'])
def analyze_batch():
    """Batch analysis for multiple sequences"""
    try:
        data = request.get_json()
        sequences = data.get('sequences', [])
        
        if not sequences:
            return jsonify({'error': 'No sequences provided'}), 400
        
        if len(sequences) > 50:  # Limit for public server
            return jsonify({'error': 'Too many sequences (maximum 50)'}), 400
        
        results = []
        for i, seq in enumerate(sequences):
            try:
                result = analyzer.predict_structure_function(seq)
                result['sequence_id'] = i + 1
                results.append(result)
            except Exception as e:
                logger.error(f"Error analyzing sequence {i+1}: {e}")
                results.append({
                    'sequence_id': i + 1,
                    'error': str(e)
                })
        
        return jsonify({
            'batch_results': results,
            'total_analyzed': len(results),
            'successful_analyses': len([r for r in results if 'error' not in r])
        })
        
    except Exception as e:
        return jsonify({'error': f'Batch analysis failed: {str(e)}'}), 500

@app.route('/api/get_educational_content')
def get_educational_content():
    """Provide educational content about structure-function relationships"""
    content = {
        'secondary_structures': {
            structure: {
                'description': analyzer.state_descriptions[structure],
                'functional_role': analyzer.functional_implications[structure]
            }
            for structure in analyzer.dssp_states
        },
        'learning_modules': [
            {
                'title': 'Chou-Fasman Method',
                'content': 'Our prediction uses the classic Chou-Fasman method with amino acid propensities.',
                'accuracy': 'Approximately 60-65% accuracy for secondary structure prediction'
            },
            {
                'title': 'Structure-Function Relationships',
                'content': 'Protein function is directly related to its 3D structure, which emerges from secondary structure.',
                'examples': ['Enzyme active sites', 'DNA-binding domains', 'Membrane channels']
            }
        ]
    }
    return jsonify(content)

if __name__ == '__main__':
    print("🧬 Simplified Protein Structure-Function Analysis")
    print("🎯 Using Chou-Fasman method for secondary structure prediction")
    print("🌐 Web interface: http://localhost:5000")
    
    app.run(debug=True, host='0.0.0.0', port=5000)