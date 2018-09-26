#ifndef PYGMO_DOCSTRINGS_HPP
#define PYGMO_DOCSTRINGS_HPP

#include <string>

namespace dcgpy
{
// kernel set
std::string kernel_init_doc(const std::string &);
std::string kernel_set_init_doc(const std::string &);
std::string kernel_set_push_back_str_doc();
std::string kernel_set_push_back_ker_doc(const std::string &);

// expression
std::string expression_init_doc(const std::string &);
std::string expression_set_doc();
std::string expression_set_f_gene_doc();
std::string expression_mutate_doc();
std::string expression_loss_doc();

// expression_weighted
std::string expression_weighted_set_weight_doc();
std::string expression_weighted_set_weights_doc();
std::string expression_weighted_get_weight_doc();

// expression_ann
std::string expression_ann_set_weight_doc();
std::string expression_ann_get_weight_doc();
std::string expression_ann_randomise_weights_doc();
std::string expression_ann_set_bias_doc();
std::string expression_ann_set_biases_doc();
std::string expression_ann_get_bias_doc();
std::string expression_ann_randomise_biases_doc();
std::string expression_ann_set_output_f_doc();
std::string expression_ann_n_active_weights_doc();
std::string expression_ann_sgd_doc();

} // namespace dcgpy

#endif
